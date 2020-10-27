# Web scraper for the iGEM wiki study

# Introduction
This tool was created for a collaboration during the iGEM 2020 contest between the [Team Manchester](2020.igem.org/Team:Manchester) and the [Team Heidelberg](2020.igem.org/Team:Heidelberg). The purpose of the web scraper is to collect data from the iGEM wiki to make it easier for future team to write and design their own wiki.

See the results [here](2020.igem.org/Team:Manchester/Wiki_Study)

## Usage

1. Install the required packages.
   The program is written in a Jupyter notebook. To open it, you need to install the Jupyter package.
   The used packages in the notebook are:
   * beautifulsoup4
   * pandas
   * numpy
   * requests
   * re

2. Run the webscraper.ipynb file.
   The file will automatically collect data from the in **return_linklist()** linked webpages. The variable **year** stores the years that will be analyzed. The **winner_(year)** list contains the winner (and for more data also the nominees) of the particular year. The for loops apply the **make_winner()** function to all entries in the list. By that, the entries in the **Winner** column are switched from False to True.

3. Use the master_dataframe_extended_winner.csv file to make your own statistical analysis.
   All statistic programs can import .csv files. Examples are: Excel, R (Studio), GraphPad Prism.
