Below are descriptions results/output files from Lester et al. 2018 Nat. Comm. deposited in Dryad as .csv files. When considering file dimensions, note that there are 1061 developable cells and that 279936 unique MSP solutions were derived.

FID1061cells.csv
1061x3 matrix of the identification number, latitude, and longitude to each of the 1061 developable cells. Latitude and longitude are in decimal degrees, using the WGS84 geodetic datum, and specify the center point for each cell.

Policy_i_a.csv
1061x279936 matrix of the 279936 MSP solutions, where 1=No Development; 2=Mussel, 3=Finfish and 4=Kelp in each of the 1061 cells

Y_NPV_wrt_MSP.csv
279936x1 vector of the dynamic response of the halibut sector to the MSP solutions

EFPayoff_a_B_wrt_DM.csv
EFPayoff_a_D_wrt_DM.csv
EFPayoff_a_F_wrt_DM.csv
EFPayoff_a_H_wrt_DM.csv
EFPayoff_a_K_wrt_DM.csv
EFPayoff_a_M_wrt_DM.csv
EFPayoff_a_V_wrt_DM.csv
Seven 1x279936 vectors, one for each sector, of the scaled (0-1), domain-wide sector value in relation to each MSP solution. The sector is indicated by the capital letter between the underscores in the file name, which represents the first letter of each sector (Mussel, Finfish, Kelp, Halibut, Viewshed, Benthos, Disease)

EFPayoff_a_ALL_wrt_DM_f01.csv 
1x279936 vector of the MSP solutions, where 0 indicates that the solution is one of the 450 filtered plans

Policy_i_a_ALL_wrt_DM_f01.csv
1061x450 matrix indicating the policy of each filtered plan, where 1=No Development; 2=Mussel, 3=Finfish and 4=Kelp

EFPayoff_a_X_wrt_DM_filter.csv 
7x450 matrix indicating the scaled (0-1), domain-wide sector values of each sector in relation to each of the filtered MSP solutions

EFPayoff_a_X_wrt_DM_Seed_bc_set.csv
7x5 matrix of the scaled (0-1), domain-wide sector values of each sector for each of the 5 seed plans