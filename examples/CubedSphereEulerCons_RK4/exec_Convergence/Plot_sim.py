import matplotlib.pyplot as plt
import pandas as pd
from numpy import sqrt 
from datetime import timedelta, datetime, date

def is_leap(Year):  
  # Checking if the given year is leap year  
    if((Year % 400 == 0) or  
     (Year % 100 != 0) and  
     (Year % 4 == 0)): 
        return True
    else:
        return False

def convert_partial_year(number):

    year = int(number)
    d = timedelta(days=(number - year)*(365 + is_leap(year)))
    day_one = datetime(year,1,1)
    date = d + day_one
    return date


def main():
    df2 = pd.read_csv('probed_data_no_interp.dat',delim_whitespace=True,skiprows=1,header=None)
    for i in range(len(df2.iloc[:,0])):
        df2.loc[i,"date"] = convert_partial_year(df2.iloc[i,0])
    df2[9] =  df2[9]/10  #microgauss to nT
    df2[10] = df2[10]/10  #microgauss to nT
    df2[11] = df2[11]/10  #microgauss to nT
    
    left = datetime(2020, 2, 6,0,0,0)
    right = datetime(2020, 2, 29,0,0,0)

    fig, ax = plt.subplots(5,1)
    ax[0].plot(df2.loc[:,"date"],df2.iloc[:,9])
    ax[0].set_ylabel('B_R (nT)')
    ax[0].get_xaxis().set_visible(False)
    ax[0].set_xlim([left, right])
    # ax[0].set_ylim([-5, 5])
    ax[1].plot(df2.loc[:,"date"],df2.iloc[:,10])
    ax[1].set_ylabel('B_T (nT)')
    ax[1].get_xaxis().set_visible(False)    
    ax[1].set_xlim([left, right])
    # ax[1].set_ylim([-5, 5])
    ax[2].plot(df2.loc[:,"date"],df2.iloc[:,11])
    ax[2].set_ylabel('B_N (nT)')
    ax[2].get_xaxis().set_visible(False)    
    ax[2].set_xlim([left, right])
    # ax[2].set_ylim([-5, 5])
    ax[3].plot(df2.loc[:,"date"],df2.iloc[:,4])
    ax[3].set_ylabel('density')
    ax[3].get_xaxis().set_visible(False)    
    ax[3].set_xlim([left, right])
    ax[4].plot(df2.loc[:,"date"],df2.iloc[:,5], label="Sim")
    plt.legend(loc="upper right")
    ax[4].set_ylabel('V_R')
    ax[4].set_xlabel('Time')
    ax[4].set_xlim([left, right])
    fig.autofmt_xdate()
    plt.show()
    # plt.savefig("probe_plot.png",dpi=300)

if __name__ == "__main__":
    main()     