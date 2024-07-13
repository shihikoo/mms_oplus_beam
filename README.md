# mms_oplus_beam

1. What is the code
   
	mms_oplus_beam is a IDL program that read data from MMS preprocessed data, identify streaming population with strong flux, narrow energy range and pointing direction. It also writes down the data into csv file, store them into tplot and make tplots.

	There are some parts of the code that are energetic that can be used for any data, which I also used for HPCA/HOPE data.

	The mms_oplus_beam, however, is a wrapper that has many MMS/HPCA specific code. Especially, the data reading part is highly dependent on Chris' data reading and process code.

	In general, if the name of the program has "mms" in it, then it has the MMS specified content. Either wise, I consider it generic.

3. Where to find the code:

	mms_oplus_beam code can be found at my home directory on our UNH machines (atlas etc.)	jliao/mms_o/*

	It can also be found and downloaded at Github. ( Open source with cc license :D )	https://github.com/shihikoo/mms_oplus_beam

	It uses the core libraries I wrote at 	jliao/cccat_user/mylib       	or https://github.com/shihikoo/mylib

5. How to call and adapt it to a different satellite.
   
	The main routine is "find_o_beam_mms.pro"

	"plot_o_beam_day_mms.pro" calls the "find_o_beam_mms.pro" for each day within the given time period.

	"find_o_beam_mms.pro" contains the main steps: read data, identify possible missing data and data errors. identify the regions, load all kinds of MMS/HPCA data, identify the beam, identify the dispersion, other calculations that I need, and then store data in csv, tplot, and make plots. The name of the function says what is for.

	One thing I want to mention is that, in the routine, we use a .tplot file for passing data. As a result, the name of the tplot file is very important to ensure the passed data are correct. As the naming convention is different for different satellites, I have a special routine called load_tplot_names_mms.pro that loads the tplot names of the tplot product. When trying to use this program for a different satellite, it is crucial to complete the load_tplot_names_xxx.pro for your specific satellite to ensure the code runs correctly.

	There are also energy steps and pitch angle steps as global environment variables, which should be changed for your own satellite. 
