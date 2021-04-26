# Aviator

### [https://ccb-compute2.cs.uni-saarland.de/aviator](https://ccb-compute2.cs.uni-saarland.de/aviator)

Aviator is a web-server monitoring the availability of other published web-servers.
It allows researchers to monitor their own tools or to asses if a tool they would like to
access is temporarily or permanently offline.

Aviator is composed of two modules:

### - [Tool List](https://ccb-compute2.cs.uni-saarland.de/aviator/tools): web-servers collected automatically from literature
### - [Aviator-enabled](https://ccb-compute2.cs.uni-saarland.de/aviator/aviator-enabled): web-servers manually added by their authors

The web-server URL or an API endpoint provided by the authors are queried twice per day. In addition to providing an availability overview we provide the possibility for authors to be notified if their webserver is offline for an unexpected period of time. 

To add your published web-server to Aviator a simple [API endpoint](https://ccb-compute2.cs.uni-saarland.de/aviator/aviator-enable) and a [registration](https://ccb-compute2.cs.uni-saarland.de/aviator/register) is needed. 

## License

[MIT © CCB-SB](../LICENSE)
