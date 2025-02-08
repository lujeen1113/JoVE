# JoVE
1. . ocean/bin/activate
2. open the python file BeautifiedGWO8.py 
	SpeedSlideWindow(2,2,23,48)
	1st Parameter is opt 2:Cplex;1:Dwave (you need to set opt=1 and fill in the your Dwave API Key: Line 405)
	2nd Parameter is addon
	3rd Parameter is HnDCount
	4th Parameter is LnDCount
	(refer to the spreadsheet)
3. Line 620 : name = '/home/lupae/Documents/Imperial_Work/IC_Campus/Documents/Proj2-JoVE/Python_codes/posdata_folder/posdata' + str(i) + '.txt'
   changed to your local directory
4. command line : python BeautifiedGWO8.py
