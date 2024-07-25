GaiaAlertsPy
=======================================

<p align="center">
  <a href='hthttps://readthedocs.org/projects/gaiaalertspy/?badge=latest'>
    <img src='https://readthedocs.org/projects/gaiaalertspy/badge/?version=latest' alt='Documentation Status' />
</a>
</p>


Welcome to GaiaAlertsPy a python package for data-mining photometric alerts from the [ESA Gaia Photometric Science Alerts](https://gsaweb.ast.cam.ac.uk/alerts). 

<p align="center">
  <img src="logo/Gapy_logo.png" width="400" height="400"/>
</p>



Installation 
==============
You can install GaiaAlertsPy using directly pip:
```shell
pip install GaiaAlertsPy
```
or install the source directly from the repository:
```shell
git clone https://github.com/AndyTza/GaiaAlertsPy.git
cd GaiaAlertsPy
python setup.py install
```


Quick-start Tutorial
=======================================

Start by querying the epochal photometric alerts for Gaia19asz. 

```python
import GaiaAlertsPy as gaap

target_id = "Gaia19asz"
alert_lc = gaap.GaiaAlert(target_id).query_lightcurve_alert()
```

Acknowledgments
=======================================
This repository was inspired by the use of the original repository [GaiaAlerts](https://github.com/davidwhogg/GaiaAlerts) (Hogg & Sipőcz) and extending its application to the time-domain community. If you use any resources or tools from this project, please cite us and the code used therein. 
