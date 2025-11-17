# -*- coding: utf-8 -*-
# Example exports of biosteam flowsheets into a standardized JSON format.
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/sff/blob/main/LICENSE
# for license details.


#%% Example 1: sugarcane-to-ethanol

from .export_codes import export_biosteam_flowsheet_sff

from biorefineries import sugarcane as sc
sc.load()
sys = sc.create_sugarcane_to_ethanol_system()
sys.simulate()
sys.diagram('cluster')

export_biosteam_flowsheet_sff(sys, "sc_ethanol_flowsheet.json")

