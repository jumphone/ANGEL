#15110700005@fudan.edu.cn
#Author: Feng ZHANG
#Python2.7
#Requirment: numpy

help_doc='''
##########################################################

             #      ##    #   #####   #####   #
           # #     # #   #   #       #       #
         #   #    #  #  #   #  ##   #####   #
       # # # #   #   # #   #   #   #       #
     #       #  #    ##   #####   #####   #####


A Network-based differential Gene Expression anaLysis tool

Author: Feng ZHANG, August 3rd 2017

Email: 15110700005@fudan.edu.cn

##########################################################

Usage:

$1: Network file (tsv, col1: gene, col2: gene)

$2: Ranked gene list

$3: Output path

##########################################################

'''

import sys,time
import numpy.ma as ma
import numpy as np





#################################################################################
#Script derived from 'stats' and 'statsmodels'

def ks_twosamp_greater(data1,data2):
    (data1, data2) = (ma.asarray(data1), ma.asarray(data2))
    (n1, n2) = (data1.count(), data2.count())
    n = (n1*n2/float(n1+n2))
    mix = ma.concatenate((data1.compressed(), data2.compressed()))
    mixsort = mix.argsort(kind='mergesort')
    csum = np.where(mixsort < n1, 1./n1, -1./n2).cumsum()
    if len(np.unique(mix)) < (n1+n2):
        csum = csum[np.r_[np.diff(mix[mixsort]).nonzero()[0],-1]]
    d = csum.max()
    prob = np.exp(-2*n*d**2)
    return prob


def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)


def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
    '''pvalue correction for false discovery rate

    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.

    Parameters
    ----------
    pvals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate
    method : {'indep', 'negcorr')

    Returns
    -------
    rejected : array, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : array
        pvalues adjusted for multiple hypothesis testing to limit FDR

    Notes
    -----

    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to alpha * m/m_0 where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).

    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
    fdr_by.

    '''
    pvals = np.asarray(pvals)

    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
##    elif method in ['n', 'negcorr']:
##        cm = np.sum(np.arange(len(pvals)))
##        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        reject_ = np.empty_like(reject)
        reject_[pvals_sortind] = reject
#        return reject_, pvals_corrected_
        return pvals_corrected_
    else:
#        return reject, pvals_corrected
        return pvals_corrected

#Test:
#b = fdrcorrection([0.002,0.02], alpha=0.05, method='indep', is_sorted=False)
#print list(b)

########################################################################################

print help_doc

try:
    NETWORK_DIR=sys.argv[1]
    INPUT_DIR=sys.argv[2]
    OUTPUT_DIR=sys.argv[3]
except Exception as e:
    exit()


print 'Start !!!'
t1=time.time()

fa=open(NETWORK_DIR)#open('/home/zhangfeng/disk/project/ANGEL/workspace/10090.protein.links.v10.5.txt.symbol.High')
GENE=set()
STRING={}
for line in fa:
    seq=line.rstrip().split('\t')
    p1=seq[0]
    p2=seq[1]
    GENE.add(p1)
    GENE.add(p2)
    if p1 in STRING:
        STRING[p1].add(p2)
    else:
        STRING[p1]=set([p2])
    if p2 in STRING:
        STRING[p2].add(p1)
    else:
        STRING[p2]=set([p1])
fa.close()



fi=open(INPUT_DIR)
genelist=[]
for line in fi:
    seq=line.rstrip().split('\t')
    gene=seq[0]
    genelist.append(gene)
      
     


GENENUM=len(genelist)
GENESCORE={}
score_all=[]
old=set()
rank=1
while rank <= GENENUM:
    gene = genelist[rank-1]
    score = abs((GENENUM-rank+1) -rank)
    if score not in old:
        score_all.append(score)
        old.add(score)
    if gene not in GENESCORE:
        GENESCORE[gene] = score
    else:
        GENESCORE[gene] =  max([score,GENESCORE[gene]])
    rank+=1

score_all.sort(reverse=True)
NUM=float(len(score_all))
RANK={}
for gene in GENESCORE:
    RANK[gene]=round(score_all.index(GENESCORE[gene])/NUM,2)



####################################
#Adjust P-value
pvalue_gene=[]
pvalue_list=[]
OUTPUT={}

i=1
for gene in genelist:
    print 'KS-test: ',i,'/',GENENUM,'\r',
    i+=1
    score_list = []
    ooo=[]
    if gene in STRING:
        for partner in STRING[gene]:
            try:
                sss=GENESCORE[partner]
                ooo.append( partner +':'+str(RANK[partner]) )
                score_list.append(sss)
            except Exception as e:
                pass
        try:
            pvalue = float(ks_twosamp_greater(score_all, score_list))
            track=';'.join(ooo)
            OUTPUT[gene]=[track]
            pvalue_list.append(pvalue)
            pvalue_gene.append(gene)
        except Exception as e:
            pass


fdr_list=list(fdrcorrection(pvalue_list))
i=0
while i<len(pvalue_list):
    gene=pvalue_gene[i]
    OUTPUT[gene].append(pvalue_list[i])
    OUTPUT[gene].append(fdr_list[i])
    i+=1
###############################################


fo=open(OUTPUT_DIR,'w')
fosig=open(OUTPUT_DIR+'.sig','w')
SIG=set()
siglist=[]
rank=1
for gene in genelist:
    if gene in OUTPUT:
        seq=OUTPUT[gene]
        track=seq[0]
        pvalue=seq[1]
        fdr=seq[2]
        fo.write(str(rank)+'\t'+gene+'\t'+str(track)+'\t'+str(pvalue)+'\t'+str(fdr)+'\n')
        if fdr<0.05:
            SIG.add(gene)
            siglist.append(gene)
            fosig.write(str(rank)+'\t'+gene+'\t'+str(pvalue)+'\t'+str(fdr)+'\n')
    else:
        fo.write(str(rank)+'\t'+gene+'\t'+'NONE'+'\t'+'NOTEST'+'\t'+'NOTEST'+'\n')
    rank+=1
fo.close()
fosig.close()


##########################################

#Drawing HIST


NUM=[]
half_window=int(len(genelist)/200)
step=2*half_window+1
i=half_window
while i < len(genelist):
    num=0
    for gene in genelist[i-half_window:i+half_window+1]:
        if gene in SIG:
           num+=1
    NUM.append([i+1,num])
    i+=step




fonum=open(OUTPUT_DIR+'.hist','w')
for one in NUM:
    num=one[1]
    fonum.write(str(one[0])+'\t'+str(num)+'\n')
fonum.close()


if len(SIG)>0:


    Rscript='''

    ##########################################################

                 #      ##    #   #####   #####   #
               # #     # #   #   #       #       #
             #   #    #  #  #   #  ##   #####   #
           # # # #   #   # #   #   #   #       #
         #       #  #    ##   #####   #####   #####


    # A Network-based differential Gene Expression anaLysis tool

    # Author: Feng ZHANG, August 3rd 2017

    # Email: 15110700005@fudan.edu.cn

    ##########################################################
    

    input="'''+OUTPUT_DIR+'.hist'+'''"
    #input_total="'''+OUTPUT_DIR+'''"
    a=read.table(input)
    #a_total=read.table(input_total)
    total='''+str(len(genelist))+'''
    sig='''+str(len(siglist)) +'''
    window='''+str(half_window*2+1) +'''
    fit=lm(a[,2] ~ I(a[,1]^2)+a[,1])
    cc=round(fit$coefficients[3]/(-2*fit$coefficients[2]),0)    
    #kmeans_input= dist(cbind(predict(fit), abs(a[,1]-cc)))
    #dis_to_cc=abs(a[,1]-cc)
    #kmeans_input= dist(cbind(scale(a[,2]), scale(dis_to_cc)))
    #kmeans_out= kmeans(kmeans_input,2)
    #cluster1=which(kmeans_out$cluster==1)
    #cluster2=which(kmeans_out$cluster==2)

    #dis_to_cc_1=dis_to_cc[cluster1]
    #dis_to_cc_2=dis_to_cc[cluster2]
    #if(mean(dis_to_cc_1) > mean(dis_to_cc_2)){
    #    tmp=cluster1
    #    cluster1=cluster2
    #    cluster2=tmp
    #    dis_to_cc_1=dis_to_cc[cluster1]
    #    dis_to_cc_2=dis_to_cc[cluster2]
    #    }
    #cutoff = min(dis_to_cc_2)
    #ll= cc-cutoff
    #rr= cc+cutoff
    pdf(paste0(input,'.pdf'))
    plot(a[,1],a[,2],pch=20,type='b',col='royalblue',xlab='Rank of center within a given window',ylab='Number of significant genes within a given window', xlim=c(0,total),ylim=c(0,max(a[,2])) ,main='ANGEL PLOT')
    par(new=T)
    plot(a[cluster2,1],a[cluster2,2],pch=1,type='p',col='red',xlab='',ylab='', xlim=c(0,total),ylim=c(0,max(a[,2])) ,main='')
    legend_1=paste0( 'Total:',as.character(total) , '   Sig:',as.character(sig),'   Window:', as.character(window))
    Pleft = round(cc/total,2) *100
    Pright = 100-Pleft
    par(new=T)
    plot(a[,1],predict(fit),type='l',col='red',lwd=2,xlim=c(0,total),ylim=c(0,max(a[,2])),ylab='',xlab='')
    lines(c(cc,cc),c(0,max(a[,2])/2.0),col='red',lwd=2)
    lines(c(ll,ll),c(0,max(a[,2])*2/3.0),col='red',lwd=1,lty=2)
    lines(c(rr,rr),c(0,max(a[,2])*2/3.0),col='red',lwd=1,lty=2)
    legend('top',legend=c(legend_1))
    text(as.integer(cc), max(a[,2])*10/20, pos=3, labels=paste0('Center:',as.character(cc)),col='red')
    #text(as.integer(ll), max(a[,2])*20/30, pos=3, labels=paste0('L-cut:',as.character(ll)),col='red')
    #text(as.integer(rr), max(a[,2])*20/30, pos=3, labels=paste0('R-cut:',as.character(rr)),col='red')
    text(as.integer(cc), max(a[,2])*9/20, pos=2, labels=paste0('Left:',as.character(Pleft),'%'),col='red')
    text(as.integer(cc), max(a[,2])*9/20, pos=4, labels=paste0('Right:',as.character(Pright),'%'),col='red')
    dev.off()
    #left_gene=a_total[which(a_total[,1]<=ll),c(1,2)]
    #right_gene=a_total[which(a_total[,1]>=rr),c(1,2)]
    #write.table(file=paste0(input_total,'.left_gene.tsv'),left_gene,sep='\t',quote=F,row.names=F,col.names=F)
    #write.table(file=paste0(input_total,'.right_gene.tsv'),right_gene,sep='\t',quote=F,row.names=F,col.names=F)


    '''
    fr=open(OUTPUT_DIR+'.R','w')
    fr.write(Rscript)
    fr.close()
    import subprocess
    subprocess.Popen('Rscript '+OUTPUT_DIR+'.R > '+OUTPUT_DIR+'.tmp',shell=True).wait()
    subprocess.Popen('rm '+OUTPUT_DIR+'.tmp',shell=True).wait()
else:
    print ''
    print 'No sig genes !!!'

###################################


t2=time.time()
ttt=round(t2-t1,2)
print ''
print 'Finished !!!'
print ''
print 'Time: '+str(ttt)+'s'
print ''
print '##########################################################'
print ''


