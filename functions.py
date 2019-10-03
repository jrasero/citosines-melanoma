import pandas as pd
import numpy as np
import matplotlib.pylab as plt

def bootci(data, stat=np.median, nboot=10000, replacement=True, alpha=0.05, method='pi',
           keepboot=False):
    """
        Compute the (1-alpha) confidence interval of a statistic (i.e.: mean, median, etc)
        of the data using bootstrap resampling.
                   
    """

    data = np.ravel(data)
    
    idx = np.random.randint(data.size, size=(nboot, data.size))

    # calculate the statistics for each bootstrap sample and sort them
    sorted_stat = np.sort(stat(data[idx], axis=1))

    # Percentile Interval method (for the moment the only one available) 
    
    ci = (sorted_stat[np.round(nboot*alpha/2).astype(int)], 
          sorted_stat[np.round(nboot*(1-alpha/2)).astype(int)])

    return ci

# Cliff's Delta from mann-whitney statistic
def delta(u, N_1, N_2): 
    return (1.0 - 2.0*u/(N_1*N_2))

# Function for the radar plots showing the performances
def plot_radar_cases(list_res_cv):
    
    radar_df = [pd.DataFrame(model) for model in list_res_cv]
    
    radar_df_mean = pd.concat([pd.DataFrame(df.loc[:,['test_bac', 'test_rec', 'test_prec','test_roc', 'test_f1w']].mean(axis=0)) 
           for df in radar_df], axis=1)
    radar_df_mean.columns = ['Breslow Thickness', 'Serum', 'Breslow Thickness + Serum']


    radar_df_std = pd.concat([pd.DataFrame(df.loc[:,['test_bac', 'test_rec', 'test_prec','test_roc', 'test_f1w']].std(axis=0)) 
           for df in radar_df], axis=1)
    radar_df_mean.columns = ['Breslow Thickness', 'Serum', 'Breslow Thickness + Serum']
    
    # ------- PART 1: Create background
    
    # number of variable
    categories = ['Balanced accuracy', 'Recall', 'Precision', 'ROC_AUC', 'F1 Weighted']
    N = radar_df_mean.shape[0]
     
    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
     
    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)
     
    # If you want the first axis to be on top:
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
     
    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, size=15, fontweight='bold')
     
    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0.6,0.7,0.8,0.9, 1.0], ["60","70","80","90","100"], color="black", size=8)
    plt.ylim(0,1)
    
    # ------- PART 2: Add plots
     
    # Plot each individual = each line of the data
    # I don't do a loop, because plotting more than 3 groups makes the chart unreadable
     
    # Ind1
    values=radar_df_mean.iloc[:,0].values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle='dotted', label="Breslow Thickness", color='Blue')
    
    y1=(radar_df_mean.iloc[:,0].values.flatten() - radar_df_std.iloc[:,0].values.flatten()).tolist()
    y1 += y1[:1]
    
    y2=(radar_df_mean.iloc[:,0].values.flatten() + radar_df_std.iloc[:,0].values.flatten()).tolist()
    y2 += y2[:1]
    
    #ax.fill_between(x=angles, y1=y1, 
     #               y2=y2, alpha=0.2, color='blue')
    
    # Ind1
    values=radar_df_mean.iloc[:,1].values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle='dashed', label="Serum", color='Royalblue')
    
    y1=(radar_df_mean.iloc[:,1].values.flatten() - radar_df_std.iloc[:,1].values.flatten()).tolist()
    y1 += y1[:1]
    
    y2=(radar_df_mean.iloc[:,1].values.flatten() + radar_df_std.iloc[:,1].values.flatten()).tolist()
    y2 += y2[:1]
    
    values=radar_df_mean.iloc[:,2].values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle='solid', label="Breslow Thickness + Serum", color='Midnightblue')
    
    y1=(radar_df_mean.iloc[:,2].values.flatten() - radar_df_std.iloc[:,2].values.flatten()).tolist()
    y1 += y1[:1]
    
    y2=(radar_df_mean.iloc[:,2].values.flatten() + radar_df_std.iloc[:,2].values.flatten()).tolist()
    y2 += y2[:1]
    
    legend = plt.legend(loc = (1.05,0.1), prop = {'weight':'bold','size':15}, frameon=True, 
               fancybox=True, edgecolor='black', title=r'${\bf Variable Domain}$')
    plt.setp(legend.get_title(),fontsize=20, color='red')
    plt.title("Logistic Regression", loc="center", pad=25, size=25, fontweight='bold')
    plt.tight_layout()
    
    return ax
#plt.savefig('scores_radar.svg')
    
# Function showing if feature has been selected in the feature selection process
def participation_plot(res_cv):
    
    n_folds = float(len(res_cv['estimator']))
    participation= 100*np.array([estim.best_estimator_.named_steps['feat'].get_support() 
                         for estim in res_cv['estimator']]).sum(axis=0)/50.
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    ax1.imshow(1-np.array([estim.best_estimator_.named_steps['feat'].get_support() \
                           for estim in res_cv['estimator']]).T[np.argsort(participation)[::-1],:], aspect='auto', cmap='gray')
    
    ax1.set_xlabel("Folds", size=15, fontweight='bold')
    ax1.set_xticks([0,9,19,29,39,49])
    ax1.set_xticklabels(["1", "10", "20","30","40","50"], size=12, fontweight='bold')
    
    ax1.set_yticks(np.arange(X.shape[1]))
    ax1.set_yticklabels(np.array(['Breslow', 'GM-CSF', 'IL-4', 'IL-6', 'IL-10', 
                                   'IL-17A', r'IFN-$\gamma$', r'TGF-$\beta$', 'DCD', 'AGE', 'SEX'])[np.argsort(participation)[::-1]],
                   size=12, fontweight='bold')
    ax1.set_yticks(0.5 + np.arange(X.shape[1]), minor=True)
    ax1.set_xticks(0.5 + np.arange(50), minor=True)
    ax1.grid(which='minor')
    ax2= ax1.twinx()
    
    ax2.set_ylim([0,9])
    ax2.set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5], minor=False)
    ax2.set_yticklabels(participation[np.argsort(participation)],
                   size=12, fontweight='bold', minor=False)
    ax1.autoscale(False)
    ax2.autoscale(False)
    ax2.set_ylabel("Total participation (%)", size=15, fontweight="bold")
    plt.title("CV selection", size=25, fontweight="bold")
    plt.tight_layout()
    
    return fig

    #plt.savefig('participation_cv.svg')


