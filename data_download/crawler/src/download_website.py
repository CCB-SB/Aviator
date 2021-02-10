import argparse
import csv
import sys
import json
import gzip
import subprocess
import traceback

from datetime import datetime
from time import sleep
from os import makedirs, rename
from os.path import join, dirname, exists

import pandas as pd

from pyvirtualdisplay import Display

from selenium import webdriver
from selenium.common.exceptions import TimeoutException, WebDriverException
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities

from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input-csv", help="Tab separated file with one URL and ID per line to check.", required=True,
)
parser.add_argument(
    "--output-folder", help="Folder where to save the results, one folder per URL", required=True,
)
parser.add_argument(
    "--output-csv", help="Tab separated file containing the results of the check", required=True,
)
parser.add_argument("--proxy-server", help="If we are behind a proxy (without authentication)")


args = parser.parse_args()

urls = pd.read_csv(args.input_csv, sep="\t")


def create_driver(log, proxy, headless=False):
    chrome_options = webdriver.ChromeOptions()

    # enable browser logging
    d = DesiredCapabilities.CHROME.copy()
    d["goog:loggingPrefs"] = {"performance": "ALL"}
    d["acceptInsecureCerts"] = True

    # chrome_options.add_argument('--headless')
    chrome_options.add_argument("--window-size=1920,1080")
    chrome_options.add_argument("--enable-precise-memory-info")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-popup-blocking")

    if proxy:
        chrome_options.add_argument(f"--proxy-server={proxy}")

    if headless:
        chrome_options.add_argument("--headless")

    driver = webdriver.Chrome(
        options=chrome_options, service_args=["--verbose", f"--log-path={log}"], desired_capabilities=d,
    )

    driver.set_page_load_timeout(30)
    driver.implicitly_wait(10)
    driver.set_window_size(1920, 1080)
    return driver


# url in response is the same as browser url
def find_response1(logs, url):
    for log in logs:
        msg = json.loads(log["message"])["message"]
        method = msg.get("method", "")
        if method == "Network.responseReceived":
            resp = msg["params"]["response"]
            msg_url = resp["url"]
            if msg_url == url:
                return {
                    "ip": resp["remoteIPAddress"],
                    "status": resp["status"],
                    "security": resp["securityState"],
                    "response": resp,
                }
    return None


# fetch document via navigation mode
def find_response2(logs):
    for log in logs:
        msg = json.loads(log["message"])["message"]
        method = msg.get("method", "")
        if method == "Network.responseReceived":
            resp = msg["params"]["response"]
            if "requestHeaders" in resp:
                sec_fetch_dest = resp["requestHeaders"].get(
                    "sec-fetch-dest", resp["requestHeaders"].get("Sec-Fetch-Dest", None)
                )
                sec_fetch_mode = resp["requestHeaders"].get(
                    "sec-fetch-mode", resp["requestHeaders"].get("Sec-Fetch-Mode", None)
                )
                if (
                    sec_fetch_dest is not None
                    and sec_fetch_dest == "document"
                    and sec_fetch_mode is not None
                    and sec_fetch_mode == "navigate"
                ):
                    return {
                        "ip": resp["remoteIPAddress"],
                        "status": resp["status"],
                        "security": resp["securityState"],
                        "response": resp,
                    }
    return None


# document and response url is contained in browser url
def find_response3(logs, url):
    for log in logs:
        msg = json.loads(log["message"])["message"]
        method = msg.get("method", "")
        if method == "Network.responseReceived":
            resp = msg["params"]["response"]
            msg_url = resp["url"]
            if msg["params"]["type"] == "Document" and msg_url in url:
                return {
                    "ip": resp["remoteIPAddress"],
                    "status": resp["status"],
                    "security": resp["securityState"],
                    "response": resp,
                }
    return None


# document and first response received
def find_response4(logs):
    for log in logs:
        msg = json.loads(log["message"])["message"]
        method = msg.get("method", "")
        if method == "Network.responseReceived":
            resp = msg["params"]["response"]
            if (
                msg["params"]["type"] == "Document"
                and "frameId" in msg["params"]
                and "loaderId" in msg["params"]
                and "requestId" in msg["params"]
            ):
                return {
                    "ip": resp["remoteIPAddress"],
                    "status": resp["status"],
                    "security": resp["securityState"],
                    "response": resp,
                }
    return None


def get_transferred_bytes(logs, request_ids):
    total_bytes = 0
    for log in logs:
        msg = json.loads(log["message"])["message"]
        method = msg.get("method", "")
        if method == "Network.dataReceived":
            total_bytes += msg["params"]["encodedDataLength"]
            if "requestId" in msg["params"]:
                request_ids.add(msg["params"]["requestId"])
    return total_bytes


def get_additional_info(driver, wait):
    status = "NA"
    ip = "NA"
    security = "NA"
    request_ids = set()
    done = False
    response = None
    transferred_bytes = 0

    url = get_current_url(driver, 5)
    all_logs = []
    for _ in range(wait):
        logs = [l for l in driver.get_log("performance")]
        all_logs.extend(logs)
        resp = find_response1(logs, url)
        transferred_bytes += get_transferred_bytes(logs, request_ids)

        if resp is None:
            resp = find_response2(logs)
            if resp is None:
                resp = find_response3(logs, url)
                if resp is None:
                    resp = find_response4(logs)

        if resp is not None:
            ip = resp["ip"]
            status = resp["status"]
            security = resp["security"]
            response = resp["response"]
            # accepted but not completed, typical for shiny app servers.
            # Check again to see if we get a different status in X (wait) seconds
            if resp["status"] != 202:
                done = True

        if done:
            break
        sleep(1)

    for l in all_logs:
        l["message"] = json.loads(l["message"])

    return {
        "code": status,
        "ip": ip,
        "certificateSecurityState": security,
        "transferred_bytes": transferred_bytes,
        "num_requests": len(request_ids),
        "response": response,
        "logs": all_logs,
    }


def get_current_url(driver, retries):
    for _ in range(retries):
        try:
            return driver.current_url
        except WebDriverException as e:
            pass
    return "NA"


def get_performance_js_measures(driver):
    navigationStart = driver.execute_script("return window.performance.timing.navigationStart")
    responseStart = driver.execute_script("return window.performance.timing.responseStart")
    domComplete = driver.execute_script("return window.performance.timing.domComplete")
    used_heap_size = driver.execute_script("return window.performance.memory.usedJSHeapSize")
    total_heap_size = driver.execute_script("return window.performance.memory.totalJSHeapSize")
    heap_size_limit = driver.execute_script("return window.performance.memory.jsHeapSizeLimit")
    transferred_bytes_page = driver.execute_script(
        "var size = 0; for(let e of window.performance.getEntriesByType('navigation')) size += e.transferSize; return size;"
    )
    transferred_bytes_resource = driver.execute_script(
        "var size = 0; for(let e of window.performance.getEntriesByType('resource')) size += e.transferSize; return size;"
    )
    num_requests = driver.execute_script(
        "return window.performance.getEntriesByType('navigation').length + window.performance.getEntriesByType('resource').length"
    )

    return {
        "backend": responseStart - navigationStart,
        "frontend": domComplete - responseStart,
        "used_heap_size": used_heap_size,
        "total_heap_size": total_heap_size,
        "heap_size_limit": heap_size_limit,
        "transferred_bytes_page": transferred_bytes_page,
        "transferred_bytes_resource": transferred_bytes_resource,
        "num_requests_determined_w_js": num_requests,
    }


def get_all_info(driver, test_url, previous_url):
    ok = False
    info = {
        "code": "NA",
        "ip": "NA",
        "certificateSecurityState": "NA",
        "transferred_bytes": "NA",
        "num_requests": "NA",
        "response": "NA",
        "logs": [],
    }
    perf = {"backend": "NA", "frontend": "NA"}
    failure_message = ""
    wait_seconds_shiny = 30

    try:
        driver.get(test_url)
        # wait at least 2 seconds
        sleep(2)
        info = get_additional_info(driver, wait_seconds_shiny)
        perf = get_performance_js_measures(driver)
        ok = True
        final_url = get_current_url(driver, 5)
    except TimeoutException as e:
        if get_current_url(driver, 5) != previous_url:
            final_url = get_current_url(driver, 5)
        else:
            final_url = test_url
        failure_message = e.msg
        print(f"timeout {e.msg}")
    except WebDriverException as e:
        final_url = get_current_url(driver, 5)
        failure_message = e.msg
        print(f"webdriver exception: {e.msg}")
        _, _, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback, limit=None, file=sys.stdout)
    return {
        "info": info,
        "perf": perf,
        "ok": ok,
        "url": final_url,
        "msg": failure_message,
    }


def check_url(driver, url):
    previous_url = get_current_url(driver, 5)
    test_url = url

    if not url.startswith("https://") and not url.startswith("http://"):
        test_url = f"http://{url}"

    res = get_all_info(driver, test_url, previous_url)

    # failure for http ? -> https
    if test_url.startswith("http://") and (not res["ok"] or res["info"]["code"] != 200):
        test_url = test_url.replace("http://", "https://")
        res2 = get_all_info(driver, test_url, previous_url)

        # replace with https if no browser failure anymore
        # or both fine and smaller error code for https
        r1code = 1e10 if res["info"]["code"] == "NA" else res["info"]["code"]
        r2code = 1e10 if res2["info"]["code"] == "NA" else res2["info"]["code"]
        if (not res["ok"] and res2["ok"]) or (res["ok"] and res2["ok"] and r2code < r1code):
            res = res2

    return {
        "url": res["url"],
        "ok": "Pass" if res["ok"] else "Fail",
        "msg": res["msg"],
        **res["info"],
        **res["perf"],
    }


def get_driver(out_dir, timestamp, retry):
    for i in range(10):
        try:
            if exists(join(out_dir, f"{timestamp}.driver.log")):
                if exists(join(out_dir, f"{timestamp}.driver.failed_{i}.log")):
                    for j in range(100):
                        if not exists(join(out_dir, f"{timestamp}.driver.failed_{i}.{j}.log")):
                            rename(join(out_dir, f"{timestamp}.driver.log"), join(out_dir, f"{timestamp}.driver.failed_{i}.{j}.log"))
                            break
                else:
                    rename(join(out_dir, f"{timestamp}.driver.log"), join(out_dir, f"{timestamp}.driver.failed_{i}.log"))
            driver = create_driver(join(out_dir, f"{timestamp}.driver.log"), args.proxy_server, retry)
            break
        except WebDriverException as e:
            # could not create chrome, retry in 2 seconds
            sleep(2)
    return driver


def process(r, out_dir, timestamp, retry=False):
    with get_driver(out_dir, timestamp, retry) as driver:
        res = check_url(driver, r["URL"])
        logs = res["logs"]
        del res["logs"]

        with open(join(out_dir, f"{timestamp}.html"), "w") as f:
            f.write(driver.page_source)
        try:
            driver.save_screenshot(join(out_dir, f"{timestamp}.png"))
        except TimeoutException as e:
            if e.msg != res["msg"]:
                res["msg"] = "save_screenshot: {}\n{}".format(res["msg"], e.msg)

        with open(join(out_dir, f"{timestamp}.cookies.json"), "w") as f:
            json.dump(driver.get_cookies(), f)

        with open(join(out_dir, f"{timestamp}.info.json"), "w") as f:
            json.dump(res, f)

        with gzip.open(join(out_dir, f"{timestamp}.logs.json.gz"), "wt", encoding="ascii") as f:
            json.dump(logs, f)

        writer.writerow(
            [
                r["ID"],
                r["URL"],
                res["url"],
                res["ok"],
                res["msg"],
                res["code"],
                res["certificateSecurityState"],
                res["transferred_bytes"],
                res["ip"],
                timestamp,
            ]
        )

    if exists(join(out_dir, f"{timestamp}.driver.log")):
        subprocess.run(["gzip", join(out_dir, f"{timestamp}.driver.log")], check=True)


with Display(visible=False, size=(1920, 1080)) as display:
    makedirs(dirname(args.output_csv), exist_ok=True)
    with open(args.output_csv, "w", newline="") as f_res:
        writer = csv.writer(f_res, dialect="excel-tab")
        writer.writerow(
            [
                "ID",
                "Original URL",
                "Derived URL",
                "Status",
                "Message",
                "Code",
                "Certificate",
                "Transfer size",
                "IP Address",
                "Timestamp",
            ]
        )

        for i, r in urls.iterrows():
            out_dir = join(args.output_folder, r["ID"])
            makedirs(out_dir, exist_ok=True)
            timestamp = datetime.now().isoformat()
            # retry up to 10 times if we crash for unknown reasons
            for i in range(10):
                try:
                    retry = i > 0
                    process(r, out_dir, timestamp, retry)
                    break
                except Exception as e:
                    print(f"Unexpected exception: {e}")
                    _, _, exc_traceback = sys.exc_info()
                    traceback.print_tb(exc_traceback, limit=None, file=sys.stdout)

