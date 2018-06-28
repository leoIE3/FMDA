#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:55:25 2017

@author: leoie

Welcome to my file :) 
"""
import numpy as np
import re
import matplotlib.pyplot as plot
from scipy import stats
from math import log10
from matplotlib import ticker

def readfile(filename):   
    """Reads the input file
    """
    with open(filename,'r') as fin:
            header=fin.readline().split('\t')
            contents=fin.readlines()
            
    """adding a primary key element and then shifting the string of contants 1 
    tab forward to have the pk before all the other data 
    """ 
    header=['Primary key']+header
    for j in range(len(contents)):
        contents[j]=contents[j].split('\t') #spliting the string into a list 
        contents[j].insert(0,j+1) #inserting a primary key into the big dataset
    return contents,header

def create_dics(data,header):
    country_code={} #country code = code:country name 
    sector_code={} #sector code = code:sector name 
    parent_sector=[] #list with parent sectors
    pollutant_name=[] #list with pollutant names
    year=[] #list with years
    for line in data:
        country_code[line[1]]=line[2]
        sector_code[line[8]]=line[6] 
        parent_sector.append(line[7])
        pollutant_name.append(line[4])
        year.append(line[5])
    return country_code,sector_code,parent_sector,pollutant_name,year


def get_children(sector_code,parent_sector):
    """Filter out parent sector codes from sector list
    """
    sectors=[] #gets a list of the sector codes in the dictionary
    
    for key in sector_code:
        sectors.append(key)   
    
    filtered_sectors=list(set(sectors)-set(parent_sector)) #narrow down the sector codes
    
    to_remove=[] #remove the sector codes starting with S or -
    for i in range(len(filtered_sectors)):
        if filtered_sectors[i].startswith('S') or filtered_sectors[i]=="-":
            to_remove.append(filtered_sectors[i])
    
    for item in filtered_sectors:
        for remove in to_remove:
            if item==remove:
                filtered_sectors.remove(remove)
                
    return sectors,filtered_sectors

def get_pollution_data(data,filtered_sectors,pollutant_name,year,country_code):
    """Extract pollution data for the chosen country
    """
    for key in country_code:
        try:
            code_to_extract=input('What is the country that you wish to extract? ')
            if code_to_extract in country_code.keys():
                break
            else:
                raise (NameError)
        except NameError:
            print('That country is unknown')
            test=input('If you want to see the list type yes: ')
            if test=='yes':
                for key in country_code:
                    print (key)
            else:
                continue
    
    sector=[] #to get the indexed sector position corresponding to the emission index
    pollution=[] #the indexed emission 
    countrypollution={} #gathers all pollution data from all countries
    for key in country_code: #initializing the dictionary
        countrypollution[key]=0
        
    for line in data:
        if code_to_extract==line[1]:
            if line[8] in filtered_sectors:
                if float(line[11]) != 0.0:
                    if line[4]!='All greenhouse gases - (CO2 equivalent)':
                        pollution.append(float(line[11]))
                        sector.append(line[8])
        else:
            if line[8] in filtered_sectors:
                if float(line[11]) != 0.0:
                    if line[4]!='All greenhouse gases - (CO2 equivalent)':
                        countrypollution[line[1]]=countrypollution[line[1]]+float(line[11])        

    return pollution,code_to_extract,countrypollution,sector

def validation(data,pollution,filtered_sectors,country_code,year,cc,countrypollution,sector,sector_code):
    """Plot country totals and sanity check
    """
    countrypollution[cc]=sum(pollution) #sums the chosen country pollution
    b=[] #just temporary stuff
    a=list(countrypollution.values())
    a.sort() #sorts the total pollution values
    if cc!='EUA' and cc!='EUC': #aim here is just to remove them from the plot because they are very large compared to others
        a.remove(countrypollution['EUA'])
        a.remove(countrypollution['EUC'])
    for item in a:
        for key in countrypollution: 
            if countrypollution[key]==item:
                b.append(key)

    fig,ax=plot.subplots()
    ax.scatter(b,a)
    ax.set(xlabel='Countries', ylabel='Sum of GG (millions)',title='Comparison of country totals')
    scale_y = 1e6
    ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
    ax.yaxis.set_major_formatter(ticks_y)
    plot.xticks(rotation=90)

    for i in range (len(a)):
        xy=(b[i],a[i])
        plot.annotate(b[i],xy)
    plot.show()
    fig.savefig("Comparison of country totals.png")

    negatives=[] #list that will have the negative items
    codes=[] #which sector have negative emissions
    for i in range(len(pollution)):
        if pollution[i]<float(0):
            negatives.append(sector[i])
    for item in negatives:
        if item in sector_code.keys():
            codes.append(sector_code[item])
    print('----------------------------------------------')
    print ('SANITY CHECK' )
    print ('The negative values occur in sectors: ' )
    print(list(set(codes)))
    print('----------------------------------------------')
              
def plot_benford(iterable,country):
    """Plot leading digit distribution in a string iterable.
    """
    numbers = [float(n) for n in range(1, 10)] #just generating a list
    benford = [log10(1 + 1 / d) for d in numbers] #calculates Benford 
    plot.rcdefaults() 
    plot.plot(numbers, benford, 'ro', label = "Benford")
    pollutionLDs = list(digits(iterable)) #yields a list with leading digits of pollution dataset
    counts=plot.hist(pollutionLDs, range(1, 11), align = 'left', normed = True,rwidth = 0.7, label = "%s"%country)
    plot.title("Benford distribution VS pollution data")
    plot.xlabel("Left-most digit")
    plot.ylabel("Frequency")
    plot.xlim(0, 10)
    plot.xticks(numbers)
    plot.legend()
    plot.savefig('histogram_%s'%country.lower())
    plot.show()
    return benford,counts[0],pollutionLDs
    
def digits(iterable):
    """Yield leading digits of number-like strings in an iterable.
    """
    numexp = re.compile(r'\d+(\.\d+)?([eE]\d+)?')
    leading = set("123456789")
    for item in iterable:
        item = str(item).replace('e-', '') #just to be sure to read also scientifical notation
        for match in numexp.finditer(item):
            for digit in match.group(0):
                if digit in leading:
                    yield int(digit)
                    break

def statistic_test(benford,counts,pollutionLDs):
    """Calculates Pearson and Chi
    """
    counts_digits = []
    counts_benford = []
    for i in range(len(benford)):
        counts_digits.append(counts[i]*len(pollutionLDs))
        counts_benford.append(benford[i]*len(pollutionLDs))
    pearson = stats.pearsonr(counts_benford, counts_digits)
    print('pearson, correlation coefficient', pearson[0])
    chi = stats.chisquare(counts_digits, counts_benford)
    print('Chi-square', chi)
    return counts_digits, counts_benford, pearson, chi 

def summary_statistics(data,cc,filtered_sectors, year): 
    """Calculates the summary statistics
    """
    years=list(set(year))
    years.remove('1985-1987') #removing this entry to not get errors because cannot integer it
    years=list(map(int,years))
    emissions_summarystats = [] #just gathering all the data - the totals
    min_year = min(years)
    max_year = max(years)
    for line in data:
        if cc==line[1]:
            if line[8] in filtered_sectors:
                if float(line[11]) != 0.0:
                    if line[4]=='All greenhouse gases - (CO2 equivalent)':
                        emissions_summarystats.append(float(line[11])) 
  
    mean = np.mean(emissions_summarystats) #calculating the mean
    median = np.median(emissions_summarystats) #calculating the median
    standarddeviation = np.std(emissions_summarystats) #calculating the SD
    skewness=stats.skew(emissions_summarystats) #calculating the skewness
    min_value = min(emissions_summarystats) 
    max_value = max(emissions_summarystats)
    normaldist=np.random.normal(mean,standarddeviation,1000)
    plot.title("Distribution of pollution data for %s, mu=%d and sigma=%d"%(cc,mean,standarddeviation))
    count, bins, ignored=plot.hist(normaldist,80,normed=True,label = 'Distribution') #plot the normal distribution
    plot.plot(bins, 1/(standarddeviation * np.sqrt(2 * np.pi)) * np.exp( - (bins - mean)**2 / (2 * standarddeviation**2) ),linewidth=2, color='r',label='PDF')
    plot.legend()
    plot.xlabel("Mean")
    plot.ylabel("Probability")
    plot.savefig('distribution_%s'%cc)
    plot.show()
    return mean, median, standarddeviation, min_value,  max_value, min_year, max_year, emissions_summarystats,skewness

def output_statistics(country, min_year, max_year, mean, median, standarddeviation, min_value, max_value, chi, pearson,skewness):
    """Writes the output file
    """
    fout = open ('statistics_%s.dat'%country.lower(), 'w') #output file with the summary of the statistics
    fout.write( 'Statistics about the emission data of {} for the time series from {} until {}\n \n'.format(country, min_year, max_year))
    fout.write( 'Name \t\t\t\t Value \t\t\t Unit\n\n')
    fout.write( 'Mean \t\t\t\t {:g} \t\t Gg CO2-eq\n'.format(mean))
    fout.write(  'Median \t\t\t\t {:g} \t\t Gg CO2-eq\n'.format(median))
    fout.write('Standard deviation \t\t {:g} \t\t Gg CO2-eq\n'.format(standarddeviation))  
    fout.write('Min. emission \t\t\t {:g} \t\t Gg CO2-eq\n'.format(min_value))
    fout.write( 'Max. emission \t\t\t {:g} \t\t Gg CO2-eq\n'.format(max_value))
    fout.write('Skewnwss \t\t {:g} \n'.format(skewness))
    fout.write( 'chi square \t\t\t {:g} \n'.format(chi[0]))
    fout.write( 'p-value of chi-square\t\t\t {:g} \n'.format(chi[1]))
    fout.write( 'Pearson correlation coefficient \t {:g} \n'.format(round(pearson[0],3)))
    fout.write('Pearson p-value \t\t\t {:g} '.format(pearson[1],3))
    fout.close() 