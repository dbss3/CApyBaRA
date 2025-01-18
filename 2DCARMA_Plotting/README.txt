####################################
2DCARMA_Plotting
Routine written by Dominic Samra c. 01/2025
Adapted from python scripts provided by Diana Powell

####################################
Basic opperation

The main routine you will need is called just_plot_it.py

To run this you just simply call it from the terminal

It has a set of required arguments though, so the call signature is:

python3 just_plot_it.py <file_path> <new_or_time_averaged_or_plotting_or_replotting=1,2,3,4>  <which_plot> <desire_value>

if you don't enter the correct number of arguments the script will tell you

<file_path> = the FULL path to a file which you will use to define the paths to 
              ALL other input and output files (see more discussion below).
              example is: ./test/files_names_and_locs.txt

<new_or_time_averaged_or_plotting_or_replotting> = integer flag used to select 
               how much of the routine needs running at a given time, this can
               save a LOT of time by avoiding re-reading in the CARMA output files,
               which is slow, though still now only on the order of ~3 mins.
               The options are:
               1 = Full run, from the raw CARMA output file
               2 = Running from the npz files of the CARMA output 
                   (i.e. starting by time averaging the data)
               3 = Running from the time averaged npz files.
                   This is the option you need if you did 1 already,
                   everything worked and now you want to maked different plots
               4 = NOT IMPLEMENTED YET: Will be a way to re-plot figures so that
                   if something went wrong with the plot (i.e. want to rescale colourbar)
                   u can using just the output arrays.
                   TBH even thinking about this now, I'd probably just make this a seperate
                   routine, too many new arguments needed/ changes depending on the changes

<which_plot> = Str option that defines what set of plots you are going to make.
               Given the current setup, you will get plots for ALL groups in the
               specified group_names_and_properties.txt file (again, more details below)
               You have no choice, why because I think you should plot everything,
               aaaand I'm too lazy to program in options for that.
               BECAUSE atm it relies entirely on ordering (there are no headers in the CARMA
               output files, working on that), so MAKE SURE the ordering matches your data
               
               Options atm are:
               lon_vs_r_numb_contour = A contourf plot of particle number density per size bin
                        with longitude on the x-axis (terminators are marked by green dashed lines)
                        and particle size in microns on the y-axis AT A SPECIFIED PRESSURE.
                        For this one therefore you must pass a float of the desired PRESSURE
                        in BAR as <desire_value>
               lon_vs_p_total_numb_contour = A contourf plot of INTEGRATED particle number density
                        i.e. the total number density of particles. Longitude is on the x-axis
                        (terminators are marked by green dashed lines) and pressure on the y-axis
                        in bar. For this one you must pass a DUMMY argument for <desire_value>,
                        I usually do 0 (has to be a float compatable type).
               r_vs_p_total_numb_contour = A contourf plot of particle number density per size bin
                        BUT it is integrated over a working angle (ATM this is fixed at +-20deg)
                        with particle radius on the x-axis (terminators are marked by green dashed lines)
                        and pressure on the y-axis in bar. For this one your <desire_value> must be a
                        float/signed int corresponding to the desired central longitude.
                        Basically, this a way to plot the terminator limbs is what I made this for.

<desire_value> = Just recapping from above, this is a float which specifies specific values that
                 have to be passed to plot at for example a certain pressure or longitude.
                 IF NO VALUE IS NEEDED MAKE SURE YOU PASS A DUMMY ARGUMENT OTHERWISE THE ROUTINE WILL ERR

####################################
####################################
file_names_and_locs file explanation:

####################################
Input
There are principly four paths that need to be constructed in order to run the routine.
Each one can be defined in this file as first a *_loc, and then *_name. Where

*_loc is the path to the file (i.e. infile) from wherever the routine is being RUN,
    EXCLUDING the last /

*_name is the name of the file, which has an assumed extension of .txt
    I.E. Don't include it here.

Each string is defined by first having a row starting with an !
e.g.:
! infile_loc
*

And then directly under that line is the input string (i.e. where * is in the above example)

In the test example two files are specified: infile and longitudes. These two MUST always be specified

Two additional files can be specified, but by default they have assumed files already.
These files define the CARMA cloud chemistry that was used so that the input can be read in and labelled
correctly.

These are:
./group_names_and_properties_DEFAULT.txt
./cloud_materials_DEFAULT.txt

If nothing is specified, these two example files will be used, 
they correspond to the chemistry used in Powell and Zhang 2024.

A couple more words on these files since they're a little jank:

./cloud_materials_DEFAULT.txt I still don't fully understand,
atm those arrays aren't used by any of the plotting scripts
it just contains the values that were in the plotting script given to me by Diana.
I'm in the process of making it work. SO for your runs this may need some fudging.

./group_names_and_properties_DEFAULT.txt contains on each line:
    the group_name
    the mass you want to assume for it for mass conversion (Not used ATM)
    the colourmap you want the contours to be in
    the elements list for the group (this entry also doesn't do anything ATM, work in progress)

These are defined by the header which is the line starting with ! DONT MESS WITH THIS LINE.
All other lines starting with # and also anything after a # on a line will be ignored.

To specify your own locations for these simply define:
cloud_properties_file_loc and cloud_properties_file_name
and/or
cloud_materials_file_loc and cloud_materials_file_name

otherwise the default value will be used.

####################################
Output
OK, last two things. For output you can specify a location to output the npz files and the plots
specify:
! outfile_loc
*
This even checks if the very last subfolder of the path exists or not and makes it if it doesn'tell

Finally, by default all npz files and plots will be tagged with the infile_name.
If you want to specify a different (maybe shorter) tag, you can do this by specifying it under
! run_name
*

####################################
####################################
Phew that's everything, if you have issues contact me (Dominic Samra) on Slack
OR by email at dominicsamra.uk@gmail.com

####################################

