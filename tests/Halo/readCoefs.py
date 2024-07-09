#!/usr/bin/env python
# coding: utf-8

import pyEXP

coefs = pyEXP.coefs.Coefs.factory('outcoef.halo.run0')
data  = coefs.getAllCoefs()
print(data.shape)
print(coefs.getName())
