import numpy as np
from tempfile import mkdtemp
import os.path

def check_listtype_imports(candidate, printname):
    # Check of data
    if not candidate:
        print "No " + printname + " supplied"
    elif type(candidate) not in [np.ndarray, list, np.core.memmap]:
        print printname + " not presented in accepted type"
    elif len(np.shape(candidate)) > 1:
        print printname + " has shape - please flatten"
    elif len(candidate) < 2:
        print "Too few points in " + printname
    else:
        return len(candidate)
    return -1


class AllanI:
    def __init__(self):
        self.status = "Init"
        self.errorflag = False
        self.data_length=0

    def importData(self, timestamps=[], rate=None, data=[], errors=[]):
        self.status="Init"
        data_length = check_listtype_imports(data, "data")
        if data_length == -1:
            return
        self.status = "Data loaded"
        self.data = data
        # Read in timestamps next
        timestamps_length = check_listtype_imports(timestamps, "timestamps")
        if timestamps_length != -1 and timestamps_length != data_length:
            print "No. Timestamps does not match No. of data points"
            return
        elif timestamps_length != -1 and timestamps_length == data_length:
            self.timestamps = timestamps
            self.status = "Timestamped Data loaded"
        else:
            if rate == None:
                print "Neither timestamps nor rate provided"
                return
            else:
                try:
                    self.rate = float(abs(rate))
                except:
                    print "Rate not provided in readable form"
                    return
        self.status = "Rate Data loaded"
        self.data_length=data_length
        # Error loading if available
        errors_length = check_listtype_imports(errors, "errors")
        if errors_length == data_length:
            self.errors = errors
            self.errorflag = True

        print "Import success"

    def setup_taus(self,tau_in=None,no_of_taus=1000,tau_min=None,minimal_stride=100,max_analysis_blocksize=1000000):
        self.tau_level=0
        if self.status=="Rate Data loaded":
            if tau_in==None:
                self.tau_list = np.logspace(np.log10(1./self.rate),np.log10(0.5*(self.data_length-1)/self.rate),int(no_of_taus))
                self.tau_list = np.unique(np.array([i for i in self.tau_list * self.rate ], dtype=long))
            else:
                try:
                    self.tau_list=np.array(tau_in)
                except:
                    print "Taulist was not provided as compatible list"
                    self.tau_list = np.logspace(np.floor(np.log10(1. / self.rate)),
                                                np.ceil(np.log10(self.data_length / self.rate)), int(no_of_taus))
                self.tau_list=np.unique(np.array([i for i in self.tau_list*self.rate if (i >= 1 and i<=self.data_length/2.)],dtype=long))
            #Stride check
            self.minimal_stride=long(minimal_stride)
            self.max_analysis_blocksize=long(max_analysis_blocksize)
            if max_analysis_blocksize > self.data_length and self.minimal_stride > self.max_analysis_blocksize/2:
                self.minimal_stride = self.max_analysis_blocksize/10
                print "Stride settings to high for data length"
            else:
                self.max_analysis_blocksize=(self.max_analysis_blocksize/self.minimal_stride)*self.minimal_stride

    def singlelevelrunner(self):
        #Make some decider about the rate thing
        self.no_of_chuncks=int(np.ceil(float(self.data_length)/self.max_analysis_blocksize))
        self.no_of_used_tau=len(self.tau_list[self.tau_list <= self.max_analysis_blocksize/self.minimal_stride])
        self.level_results=[[[] for i in range(self.no_of_used_tau)] for j in range(self.no_of_chuncks)]

        for chunck in range(self.no_of_chuncks):
            for tau in range(self.no_of_used_tau):
                if self.errorflag:
                    error=self.data[chunck * self.max_analysis_blocksize:(chunck + 1) * self.max_analysis_blocksize]
                else:
                    error=None
                self.level_results[chunck][tau]=self.singleAllan(self.data[chunck*self.max_analysis_blocksize:(chunck+1)*self.max_analysis_blocksize],
                            self.tau_list[tau],error=error)
        self.level_switcher()

    def level_switcher(self):
        #Make "means from" multiple chunks
        self.level_results=np.array(self.level_results,dtype=float)
        for tau in range(self.no_of_used_tau):
            if self.errorflag:
                weights = 1. / self.level_results[:, tau, 2] ** 2
            else:
                weights = self.level_results[:, tau, 2]  # Just the number of entries
            self.global_results[self.tau_list[tau]]=[np.average(self.level_results[:,tau,1],weights=weights),
                 1./np.sqrt(np.sum(weights)),self.tau_level]
        #create new data basis for next level
        if self.data_length>self.max_analysis_blocksize:
            self.tau_level+=1
            self.tau_list=np.unique(np.array(np.array(self.tau_list[self.no_of_used_tau:])/self.minimal_stride,dtype=long))
            new_data_dir = os.path.join(mkdtemp(), 'new_allan_data_'+str(self.tau_level)+'.dat')
            new_data = np.memmap(new_data_dir, dtype='float32', mode='w+',shape=(self.data_length/self.minimal_stride))
            if self.errorflag:
                new_error_dir = os.path.join(mkdtemp(), 'new_allan_error'+str(self.tau_level)+'.dat')
                new_errors = np.memmap(new_error_dir, dtype='float32', mode='w+',shape=(self.data_length/self.minimal_stride))
            for i in range(0,(self.data_length/self.minimal_stride)*self.minimal_stride,self.minimal_stride):
                if self.errorflag:
                    weights=1./np.array(self.errors[i:i+self.minimal_stride])**2
                    new_data[i/self.minimal_stride]=np.average(np.array(self.data[i:i+self.minimal_stride]),weights=weights)
                    new_errors[i/self.minimal_stride]=1./np.sqrt(np.sum(weights))
                else:
                    new_data[i/self.minimal_stride]=np.mean(np.array(self.data[i:i+self.minimal_stride]))
            if self.tau_level>1:
                del self.data
                self.data=new_data
                if self.errorflag:
                    del self.errors
                    self.errors=new_errors
            self.singlelevelrunner()


        else:
            self.generate_results()


    def generate_results(self):
        k=1
        #Make output processing

    def oadev(self):
        #Has to be called to start the allan deviation analysis
        self.global_results={}  #Dictionary, indexed by tau/tau_0. Filled with list. 0:2*sig**2, 1: Delta 2*sig**2 2:taulevel
        self.singlelevelrunner()


    def singleAllan(self,data,tau,error=None):
        return [tau,1.,2.]


test_obj = AllanI()
test_obj.importData(rate=1., data=[i for i in xrange(1200000)])
test_obj.setup_taus(no_of_taus=7)
print test_obj.tau_list
test_obj.oadev()
print test_obj.global_results