import bs4 as bs
import pandas as pd
import os
import zipfile
import numpy as np
import logging

logging.basicConfig(format='[%(levelname)s %(asctime)s]: %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

def change_format(x):
    try:
        return ''.join(''.join(x.strip().split('<P>')).split('  '))
    except:
        return '-'
def newline_killer(x):
    try:
        return ('').join(x.split('\n'))
    except:
        return x

if os.path.isfile("xml_parts.xml"):
    logging.info("xml_parts.xml does already exist. It will not be extracted from xml_parts.xml.zip")
else:
    logging.info("xml_parts.xml does not exists. It will be extracted now.")
    with zipfile.ZipFile("xml_parts.xml.zip", 'r') as zip_ref:
        zip_ref.extractall()

if os.path.isfile("xml_parts_good.xml"):
    logging.info("xml_parts_good.xml does already exist. It will not be created.")
else:
    logging.info("xml_parts_good.xml does not exist yet. It will be created now.")
    with open('xml_parts.xml', encoding='utf-8', errors='ignore') as f:
        content = f.read().splitlines()
        f.close()
    f = open("xml_parts_good.xml", "a")
    for item in content:
        f.write(item+"\n")
    f.close()
    logging.info("xml_parts_good creation finished.")
with open("xml_parts_good.xml", "r") as xml_file:
    logging.info("xml file will be parsed. This can take a while.")
    soup = bs.BeautifulSoup(xml_file, "html.parser")
    logging.info("xml file parsed successfully.")
columnnames = []
for hit in soup.find_all("field"):
    if hit.get("field") != None:
        columnnames.append(hit.get("field"))
columnnames = set(columnnames)

dataframe = pd.DataFrame(columns=columnnames)
total = len(soup.find_all("row"))
counter = 0
for row in soup.find_all("row"):
    dataframe_temp = pd.DataFrame(index=[0], columns=columnnames)
    for columnname in columnnames:
        try:
            dataframe_temp.at[0, columnname] = row.find(
                "field", {"name": columnname}).string
        except:
            pass
    frames = [dataframe, dataframe_temp]
    dataframe = pd.concat(frames)
    counter += 1
    if counter % 1000 == 0:
        print("Current progress:",end=" ")
        print(round((counter/total)*100, 2), "%")
logging.info("Done!")
dataframe = dataframe.replace(r'^\s*$', np.nan, regex=True)
dataframe = dataframe.dropna(axis=1, how="all")
dataframe = dataframe.dropna(axis=0, how="all")
dataframe = dataframe.drop(columns=["sequence_sha1"])
dataframe.notes = dataframe.notes.apply(change_format)
dataframe = dataframe.applymap(newline_killer)
dataframe.to_csv("registry.csv",index=False)


