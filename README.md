GaiaAlertsPy
=======================================
Welcome to `GaiaAlertsPy` a small python package for data-mining photometric alerts from the [ESA Gaia Photometric Science Alerts](https://gsaweb.ast.cam.ac.uk/alerts). 


Installation 
==============
You can install GaiaAlertsPy using directly `pip`:
```shell
pip install GaiaAlertsPy
```
or install the source directly from the repository:
```shell
pip install git+ssh://git@github.com/AndyTza/GaiaAlertsPy.git
```


Quick-start Tutorial
=======================================

Start by querying the epochal photometric alerts for Gaia19asz. 

```python
import GaiaAlertsPy as gappy

target_id = "Gaia19asz"
alert_lc = gappy.GaiaAlert(target_id).query_lightcurve_alert()
```

Acknowledgments
=======================================
This repository was inspired by the use of the original repository [GaiaAlerts](https://github.com/davidwhogg/GaiaAlerts) (Hogg & Sipőcz) and extending its application to the time-domain community. If you use any resources or tools from this project, please cite us and the code used therein. 
