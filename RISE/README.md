# RISE

# Introoduction
[RISE](2020.igem.org/Team:Heidelberg/Software/RISE) is a tool developed by the iGEM [Team Heidelberg](2020.igem.org/Team:Heidelberg) 2020. All of us use the iGEM Registry to search for specific parts to use in our projects. We saw the need for a better searchability of the registry. That's why we created RISE (**R**egistry **I**ntelligent **S**earch **E**ngine), an easy-to-use tool to search the iGEM registry. Together with the actual search engine, a tool to locally create the registry as a .csv file is provided.

## License
The code is published under the MIT license.

## Usage

1. Install the required packages.
   The requires packages are:
   * BeautifulSoup4
   * Pandas
   * Numpy
   * Biopython
   * Jupyter
   * IPython Widgets
   * Qgrid


2. Create the .csv file.
   The [.xml dump](http://parts.igem.org/Registry_API#POINT-IN-TIME_DATABASE_DUMP) (20 November 2019) is provided as a zip file in the GitHub repository. When you run the csv_creator.py script, it will be unpacked if the file xml_parts.xml does **not** exist. If a new dump is available, you can download it and put it in the RISE directory. The file will be parsed and a .csv file will be created. To make the decision if you want to use RISE easier for you, we provided the first 10.000 rows in this repository for you to play with.

3. Use RISE.ipynb to search the registry.
   After you opened the Jupyter notebook, run all cells. You should now see an interactive data frame. Use the funnel symbol on the columns to filter the data frame. After you filtered, you can select a row by clicking on it, multiple separated rows by holding command (mac) or control (windows) or multiple connected rows by holding shift.
   After you selected the rows you want, click on the "Add to export data frame" button. The rows will be added to a separate data frame. You can now filter many times and add each row to the data frame until you are finished. Afterwards, you can export to a .csv file or a .fasta file by clicking one of the buttons. The file will be stored in the working directory.
