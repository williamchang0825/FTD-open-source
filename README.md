# Systems Biology Methods
## (1) Constructing candidate PPI and GRN [Preprocess]
* Run: assignment.m
* Import: GRN_name.mat
* Set: Virus name [InVirus = {A,B,…}]
* Output: PPI_name.mat、 GRN_name.mat、PPI.mat、GRN.mat
* Requirement: PPI_EXP.mat、GRN_EXP.mat (Sort them with the name constructed)
* Run: spline.m
* Import: Expression.xlsx
* Set: timepoints [default: spline([0 240 480 1440],V,linspace(0,1440,1440))]
 Ex. If your time points is 10, 15, 20, 25 you want to spline it into 1050 points
 Then set your values as spline([0 350 700 1050],V,linspace(0,1050,1050))
## (2) Executing ID
###### PPI
* Run: ppi_ID2.m
* Import: PPI_name.mat、PPI.mat、PPI_EXP.mat
* Set: start = ??? (default: id [start from the beginning])
 m = 1:nnn/k (default m = 1:nnn/2) [larger k requires shorter running period]
* Temporary output: PPI_ID_TEMP.mat [save all the variable]
* Output: PPI_PNP.mat、PPI_Edge.mat、Basal_PPI.txt
###### GRN
* Run: grn_ID.m
* Import: GRN_name.mat、GRN.mat、GRN_EXP.mat
* Output: GRN_PNP.mat、GRN_Edge.mat、Basal_GRN.txt
* ## (3) PNP
* Run: pnp_modified.m
* Import: PPI_name.mat、 GRN_name.mat、PPI_PNP.mat、GRN_PNP.mat
* Output: BB1.mat、BB2.mat、Gup.mat、Gdn.mat
* Run: rank.py
* Import: BB1.mat、BB2.mat、Gup.mat、Gdn.mat
* Output: Rank.csv
