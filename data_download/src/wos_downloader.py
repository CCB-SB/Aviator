from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options

import shutil
import time
import os
import os.path
from threading import Thread

from tqdm import tqdm

import random
import math

import pandas as pd

import sys

def sleep():
    time.sleep(random.random()+ 0.7)

num_threads = 1
results_dir = sys.argv[2]

pubs = pd.read_csv(sys.argv[1], sep='\t')
ids = pubs["PMID"].astype(str)
query = "PMID=( {} )".format(" OR ".join(ids))

os.makedirs(results_dir, exist_ok=True)

with open(os.path.join(results_dir, "query.txt"), 'w') as f:
    f.write(query)


profile = webdriver.FirefoxProfile()
profile.set_preference('browser.download.folderList', 2)
profile.set_preference('browser.download.manager.showWhenStarting', False)
profile.set_preference('browser.download.dir', os.path.abspath(results_dir))
profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/csv,text/plain')

drv = webdriver.Firefox(profile)
wait = WebDriverWait(drv, 10)
drv.implicitly_wait(10)
drv.get("https://apps.webofknowledge.com/")
elem = drv.find_element_by_link_text("Advanced Search")
elem.click()
sleep()

#query_box.send_keys(query)

sleep()
try:
    cookie_button = drv.find_element_by_xpath('//*[@id="onetrust-reject-all-handler"]')
    cookie_button.click()
except:
    print("No Cookie-Button")
sleep()
query_box = drv.find_element_by_xpath('//*[@id="advancedSearchInputArea"]')
query_box.click()
sleep()
s1 = query[:5]
for c in s1:
    query_box.send_keys(c)
sleep()
sx = query[:-1]
drv.execute_script("arguments[0].value = arguments[1];", query_box, sx)
sleep()
query_box.send_keys(")")
sleep()
search_button = drv.find_element_by_xpath('//*[@cdxanalyticscategory="wos-search-search-adv"]')
search_button.click()

sleep()

#results_link = drv.find_element_by_xpath('//*[@id="hitCount"]')
results_link = drv.find_element_by_xpath('//*[@class="brand-blue"]')
#results_link.click()

sleep()
sleep()
sleep()

#while drv.find_element_by_id("hitCount.top").text == "...":
#    time.sleep(2)

num_records = int(results_link.text.replace(',', ''))
print(num_records)

def download_records(start, end, num_records):
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2)
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    profile.set_preference('browser.download.dir', os.path.abspath(results_dir))
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/csv,text/plain')

    options = Options()

    driver = webdriver.Firefox(firefox_profile=profile, options=options)

    wait = WebDriverWait(driver, 10)

    driver.implicitly_wait(10)

    driver.get("https://apps.webofknowledge.com/")

    elem = driver.find_element_by_link_text("Advanced Search")
    elem.click()

    sleep()
    sleep()
    try:
        cookie_button = driver.find_element_by_xpath('//*[@id="onetrust-reject-all-handler"]')
        cookie_button.click()
    except:
        print("No Cookie-Button")
    sleep()
    ################################
    query_box = driver.find_element_by_xpath('//*[@id="advancedSearchInputArea"]')
    query_box.click()
    sleep()
    s1 = query[:5]
    for c in s1:
        query_box.send_keys(c)
    sleep()
    sx = query[:-1]
    driver.execute_script("arguments[0].value = arguments[1];", query_box, sx)
    sleep()
    query_box.send_keys(")")
    sleep()
    search_button = driver.find_element_by_xpath('//*[@cdxanalyticscategory="wos-search-search-adv"]')
    search_button.click()
    ################################
    #query_box = driver.find_element_by_xpath('//*[@id="value(input1)"]')
    #driver.execute_script("arguments[0].value = arguments[1];", query_box, query)
    #query_box.send_keys(query)

    #sleep()
    #time.sleep(10)

    #search_button = driver.find_element_by_xpath('//*[@id="search-button"]')
    #search_button.click()

    sleep()

    #results_link = driver.find_element_by_xpath('//*[@id="hitCount"]')
    #results_link.click()
    results_link = driver.find_element_by_xpath('//*[@class="brand-blue"]')

    sleep()

    #WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="exportTypeName"]')))
    WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.XPATH, '//*[@class="mat-icon notranslate svg-icon-18 margin-left-3--reversible mat-icon-inline mat-icon-no-color"]')))

    #while driver.find_element_by_id("hitCount.top").text == "...":
    #    sleep()

    #num_records_new = int(driver.find_element_by_id("hitCount.top").text.replace(',', ''))
    num_records_new = int(results_link.text.replace(',', ''))
    assert num_records_new == num_records

    #export_menu_button = driver.find_element_by_xpath('//*[@id="exportTypeName"]')
    export_menu_button = driver.find_element_by_xpath('//*[@class="mat-icon notranslate svg-icon-18 margin-left-3--reversible mat-icon-inline mat-icon-no-color"]')
    export_menu_button.click()
    sleep()

    #save_to_menu = driver.find_element_by_id('saveToMenu')
    #other_file_format_item = save_to_menu.find_element_by_xpath('li[4]/a')
    #other_file_format_item.click()
    save_to_menu = driver.find_element_by_id('exportToFieldTaggedButton')
    #save_to_menu = driver.find_element_by_id('exportToBibtexButton')
    save_to_menu.click()
    sleep()

    def export_records(start, end, is_first_call):
        if not is_first_call:
            #WebDriverWait(driver, 20).until(EC.invisibility_of_element_located((By.XPATH, "//div[@class='outputOptionsMsgTemplate']")))
            sleep()
            sleep()
            sleep()
            sleep()
            sleep()
            sleep()
            sleep()
            sleep()
            #export_menu = WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.XPATH, '//*[@class="mat-focus-indicator cdx-but-md mat-stroked-button mat-button-base mat-primary"]')))
            #export_menu.click()
            export_menu_button = driver.find_element_by_xpath('//*[@class="mat-icon notranslate svg-icon-18 margin-left-3--reversible mat-icon-inline mat-icon-no-color"]')
            export_menu_button.click()
            sleep()
            save_to_menu = driver.find_element_by_id('exportToFieldTaggedButton')
            save_to_menu.click()
            sleep()
        record_from_to_input_radio = wait.until(EC.element_to_be_clickable((By.ID, 'radio3')))
        record_from_to_input_radio = driver.find_element_by_xpath('//*[@for="radio3-input"]')
        record_from_to_input_radio.click()
    #    sleep()
        #record_from = driver.find_element_by_id("mat-input-2")
        record_from = driver.find_element_by_xpath('//*[@aria-label="Input starting record range"]')
        
        #record_to = driver.find_element_by_id("mat-input-3")
        record_to = driver.find_element_by_xpath('//*[@aria-label="Input ending record range. A maximum of 1000 records can be exported at one time."]')
        record_from.clear()
        record_from.send_keys(start)
    #    sleep()
        record_to.clear()
        record_to.send_keys(end)
    #    sleep() 
        #record_content_select_menu = driver.find_element_by_xpath('//*[@class="quick-output-section-content"]/span[contains(@class, "select2-container")]')
        record_content_select_menu = driver.find_element_by_xpath('//*[@class="margin-top-5 ng-star-inserted"]')
        record_content_select_menu.click()
    #    sleep()
        #full_record_option = driver.find_element_by_xpath("//*[@id='select2-bib_fields-results']/li[3]")
        full_record_option = driver.find_element_by_xpath('//*[@title="Full Record"]')
        full_record_option.click()
        sleep()
        #export_btn = driver.find_element_by_id("exportButton")
        export_btn = driver.find_element_by_xpath('//*[@class="mat-focus-indicator cdx-but-md mat-stroked-button mat-button-base mat-primary"]')
        export_btn.click()
        sleep()
        try:
            close_btn = driver.find_element_by_partial_link_text("Close")
            close_btn.click()
        except:
            pass
        shutil.move(os.path.join(results_dir, 'savedrecs.txt'), os.path.join(results_dir, 'savedrecs_{}-{}.txt'.format(start,end)))

        if is_first_call:
            # this blocks the clicking -> so simply hide it
            # waiting takes several seconds
            i = 1
            while True:
                try:
                    element = driver.find_element_by_xpath("//div[contains(@class, 'ui-dialog-titlebar') and contains(@class, 'ui-widget-header')]")
                    #driver.execute_script("arguments[0].style.visibility='hidden'", element)
                    driver.execute_script("arguments[0].remove()", element)
                    print("remove {}".format(i))
                    i += 1
                except:
                    break
            #element = driver.find_element_by_class_name("outputOptionsMsgTemplate")
            #driver.execute_script("arguments[0].style.visibility='hidden'", element)
    #        wait.until(EC.invisibility_of_element_located((By.XPATH,
    #              "//div[@class='ui-dialog-titlebar ui-widget-header ui-helper-clearfix']")))


    first_call = True
    for i in tqdm(range(start, end+1, 500)):
        print("export {}-{}".format(i, min(i+499, end),))
        
        try:
            export_records(i, min(i+499, end), i == 1)
            first_call = False
        except Exception as e:
            print(e)
            print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
            driver.quit()
            profile = webdriver.FirefoxProfile()
            profile.set_preference('browser.download.folderList', 2)
            profile.set_preference('browser.download.manager.showWhenStarting', False)
            profile.set_preference('browser.download.dir', results_dir)
            profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/csv,text/plain')
            driver = webdriver.Firefox(profile)
            wait = WebDriverWait(driver, 10)

            driver.implicitly_wait(10)

            driver.get("https://apps.webofknowledge.com/")

            elem = driver.find_element_by_link_text("Advanced Search")
            elem.click()

            sleep()
            try:
                cookie_button = driver.find_element_by_xpath('//*[@id="onetrust-reject-all-handler"]')
                cookie_button.click()
            except:
                print("No Cookie-Button")
            sleep()

            #query_box = driver.find_element_by_xpath('//*[@id="value(input1)"]')
            query_box = driver.find_element_by_xpath('//*[@id="advancedSearchInputArea"]')
            #driver.execute_script("arguments[0].value = arguments[1];", query_box, query)
            query_box.click()
            sleep()
            s1 = query[:5]
            for c in s1:
                query_box.send_keys(c)
            sleep()
            sx = query[:-1]
            driver.execute_script("arguments[0].value = arguments[1];", query_box, sx)
            sleep()
            query_box.send_keys(")")
            sleep()
            search_button = driver.find_element_by_xpath('//*[@cdxanalyticscategory="wos-search-search-adv"]')
            search_button.click()

            #search_button = driver.find_element_by_xpath('//*[@id="search-button"]')
            #search_button.click()

            sleep()
            sleep()
            sleep()

            #results_link = driver.find_element_by_xpath('//*[@id="hitCount"]')
            #results_link.click()

            #sleep()
            #while driver.find_element_by_id("hitCount.top").text == "...":
            #    sleep()
            
            results_link = driver.find_element_by_xpath('//*[@class="brand-blue"]')

            #num_records_new = int(driver.find_element_by_id("hitCount.top").text.replace(',', ''))
            num_records_new = int(results_link.text.replace(',', ''))
            assert num_records == num_records_new

            #export_menu_button = driver.find_element_by_xpath('//*[@id="exportTypeName"]')
            export_menu_button = driver.find_element_by_xpath('//*[@class="mat-icon notranslate svg-icon-18 margin-left-3--reversible mat-icon-inline mat-icon-no-color"]')
            export_menu_button.click()
            sleep()

            #save_to_menu = driver.find_element_by_id('saveToMenu')
            #other_file_format_item = save_to_menu.find_element_by_xpath('li[4]/a')
            save_to_menu = driver.find_element_by_id('exportToBibtexButton')
            #other_file_format_item.click()
            save_to_menu.click()
            sleep()
            first_call = True
    driver.quit()


records_per_thread = math.ceil(num_records/num_threads)
threads = [Thread(target=download_records, args=(start, min(start + records_per_thread-1, num_records), num_records)) for start in range(1, num_records+1, records_per_thread)]

for t in threads:
    t.start()

for t in threads:
    t.join()
