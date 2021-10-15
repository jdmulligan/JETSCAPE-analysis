# Settings in JETSCAPE
## Binning of pthat (in GeV)
[ [1,2], [2,3], [3, 4], [4, 5], [5,7], [7,9], [9,11], [11,13], [13,15], [15,17], [17,20], [20,25], [25, 30], [30, 35], [35, 40], [40, 45], [45, 50], [50, 55], [55, 60], [60, 70], [70,80], [80,90], [90, 100], [100, 110], [110, 120], [120, 130], [130, 140], [140, 150], [150, 160], [160, 170], [170, 180], [180, 190], [190, 200], [200, 210], [210, 220], [220, 230], [230, 240], [240, 250], [250, 260], [260, 270], [270, 280], [280, 290], [290, 300], [300, 350], [350, 400], [400, 450], [450, 500], [500, 550], [550, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000,1100], [1100,1200], [1200,1300], [1300,1400], [1400,1500], [1500,1600], [1600,1700], [1700,1800], [1800,1900], [1900,2000], [2000,2200], [2200,2400], [2400,2510] ]

Below is for yaml
  pt_hat_bins:
    label: 'pt_hat_bins'
    values: [1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510]

## Total Event Number for each pthat bin
120000

## Models

MATTER for vacuum (pp)

MATTER and MATTER+LBT for in-medium (PbPb)

## Parameters

alpha_S = 0.25, 0.2, 0.3

Q_switch = 1, 2, 3 GeV for MATTER+LBT

Recoil is ON


## Centralities
0-10%, 30-40%, 40-50%

# Structures of Data Files (ASCII)

## JetscapeHadronList
0:id 1:pid 2:particle status (-1 for negative particles) 3:E 4:Px 5:Py 6:Pz 7:Eta 8:Phi

## SigmaHard
0:sigma 1:error
(both are accumulated)

## SoftAnisotropy_Coefficients
0:event id 1:v2 reference 2:psi2 reference 3:v3 reference 4:psi3 reference

# Files and Folders

All the files and folders are in `AAPaperData` folder on OSIRIS (JETSCAPE-NEW)

pp(vacuum), MATTER:
`5020_PP`

PbPb, centrality 0-10%, alpha_S = 0.25, Q_switch=1.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_0-10_0R25_1R0_1  `

PbPb, centrality 0-10%, alpha_S = 0.25, Q_switch=2.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_0-10_0R25_2R0_1`  

PbPb, centrality 0-10%, alpha_S = 0.25, Q_switch=3.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_0-10_0R25_3R0_1` 

PbPb, centrality 0-10%, alpha_S = 0.2, Q_switch=2.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_0-10_0R2_2R0_1` 

PbPb, centrality 0-10%, alpha_S = 0.3, Q_switch=2.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_0-10_0R2_2R0_1` 

PbPb, centrality 0-10%, alpha_S = 0.25, MATTER, recoil ON:  
`5020_PbPb_0-10_0R25_-1R0_1` 



-
PbPb, centrality 30-40%, alpha_S = 0.25, Q_switch=1.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_30-40_0R25_1R0_1`  

PbPb, centrality 30-40%, alpha_S = 0.25, Q_switch=2.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_30-40_0R25_2R0_1`  

PbPb, centrality 30-40%, alpha_S = 0.25, Q_switch=3.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_30-40_0R25_3R0_1`  

PbPb, centrality 30-40%, alpha_S = 0.25, MATTER, recoil ON:
`5020_PbPb_30-40_0R25_-1R0_1` 

-
PbPb, centrality 40-50%, alpha_S = 0.25, Q_switch=1.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_40-50_0R25_1R0_1`  

PbPb, centrality 40-50%, alpha_S = 0.25, Q_switch=2.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_40-50_0R25_2R0_1`

PbPb, centrality 40-50%, alpha_S = 0.25, Q_switch=3.0 GeV (MATTER+LBT), recoil ON:
`5020_PbPb_40-50_0R25_3R0_1`

PbPb, centrality 40-50%, alpha_S = 0.25, MATTER, recoil ON:
`5020_PbPb_40-50_0R25_-1R0_1` 






