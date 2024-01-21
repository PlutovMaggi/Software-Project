from warnings import WarningMessage
import spkmeansmodule #prevously - mykmeanssp
import sys
import numpy as np
import pandas as pd

def kmeans():
    #recieve and check arguments
    if(len(sys.argv)!=4):
        sys.exit("Invalid Input!")

    #get arguments
    k = sys.argv[1]
    cur_goal = sys.argv[2]
    input = sys.argv[3]

    if (not assert_input(k,cur_goal,input)): #asssert input argument
        sys.exit("Invalid Input!")
        
    k = int(k)
    max_iter = 300
    epsilon = 0

    #get datapoints
    try:
        input_file = pd.read_csv(input,header=None)
    except:
        sys.exit("Invalid Input!") #file does not exists, problem reading file
    
    try:
        datapoints = input_file.to_numpy(dtype = float)
        (n,d) = datapoints.shape

        if(k<0 or k>=n): #assert k
            sys.exit("Invalid Input!")
        
        datapoints = datapoints.tolist()
        
        finalresult = spkmeansmodule.cases(n,k,d,cur_goal,datapoints) #run the c extenstion
        
        if (cur_goal == "spk"): #continue with step 6 of algorithm
            datapoints = finalresult
            d = len(finalresult[0])
            #update calculated k in the c extension
            if (k==0): 
                k = d #number of columns in first datapoint
            finalresult = run_spk(n,k,d,datapoints,max_iter,epsilon)
        
        #print final results
        print_results(finalresult)

    except SystemExit:
        sys.exit("Invalid Input!")     
    except:
        sys.exit("An Error Has Occurred")

#help function - assert input arguments
#returns false in case of an error
def assert_input(k,cur_goal,input):
    #check cur_goal values
    if(cur_goal!="spk" and cur_goal!="wam" and cur_goal!="ddg" and cur_goal!="lnorm" and cur_goal!="jacobi"):
        return False

    #check input file 
    ends = [".txt",".csv"]
    if((not k.isnumeric()) or (input[-4:] not in ends)): #asset for other parameters
        return False
    
    return True

#help function - run & print spk case
def run_spk(n,k,d,datapoints,max_iter,epsilon):
    #initialize centroids
    centroids,selections = kmeanspp(k,d,datapoints) 
    centroids = centroids.tolist()

    #run the c extenstion
    centroids = spkmeansmodule.fit(n,d,k,max_iter,epsilon,centroids,datapoints)

    if centroids == None:
        raise Exception()

    #print results
    final_selections = ','.join(str(selection) for selection in selections)
    print (final_selections)
    #return final result
    return centroids;

#kmeans pp from hw2
def kmeanspp(k,d,datapoints):
    selections = []
    np.random.seed(0)
    N = len(datapoints)
    i = 1
    m = np.random.choice(N) #choose random index
    selections.append(m)
    centroids = np.zeros((k,d))
    centroids[0] = np.copy(datapoints[m])
    while i<k:
        D = np.zeros(N)
        P = np.zeros(N)
        for l in range(N):
            mindelta = float('inf')
            for j in range(i):
                delta = np.sum((datapoints[l]-centroids[j])**2)
                if(delta<mindelta ):#new min
                    mindelta = delta
            D[l] = mindelta
        
        sumOfD = np.sum(D) 
        if sumOfD!=0:   
            for l in range(N):    
                P[l] = D[l]/sumOfD
                
        m = np.random.choice(N,p=P)
        selections.append(m)
        centroids[i] = np.copy(datapoints[m])
        i+=1
    return centroids,selections

#print final result matrix
def print_results(finalresult):
    final_matrix = ""
    for i in range(len(finalresult)):
        for j in range(len(finalresult[0])):
            final_matrix += "%.4f" % finalresult[i][j]
            if(j != (len(finalresult[0])-1)):
                final_matrix +=","
        if (i != (len(finalresult)-1)):
            final_matrix += "\n"
    print(final_matrix)

if __name__ == '__main__':
    kmeans()
 