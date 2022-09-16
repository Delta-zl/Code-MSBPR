# MSBPR
coda for MSBPR(MSBPR: A Multi-pairwise Preference and Similarity Based Bayesian Personalized Ranking Method for Recommendation)

## Datasets Information
30Music:   -n 4586 -m 4989 -https://recsys.deib.polimi.it/?s=30music&submit=Search  
Tmall:	   -n 4516 -m 5000 -https://tianchi.aliyun.com/datalab/dataSet.html?    
rsc15:     -n 4926 -m 4997 -http://2015.recsyschallenge.com/challenge.html  
aotm:      -n 4885 -m 5000 -https://bmcfee.github.io/data/aotm2011.html  
Last.fm:   -n 1640 -m 12440 -https://grouplens.org/datasets/hetrec-2011/  
Yelp:      -n 4978 -m 5000  -https://www.yelp.com/dataset  

## Configurations for these datasets
### 30music-5k5k
fnTrainData = "./data/30music/30music-5k5k-train";  
fnTestData = "./data/30music/30music-5k5k-test";  
fnValidData = "./data/30music/30music-5k5k-test";  
fnSimilarity_i = "./similarity_matrix/30music/similarity_i.txt";  
fnSimilarity_u = "./similarity_matrix/30music/similarity_u.txt";  
fnSimilarity_i_u = "./similarity_matrix/30music/similarity_i_u.txt";  
n = 4586, m = 4989 ;  
alpha_u = alpha_v = beta_v = 0.1;  
gamma = 0.01;  
threshold_p = 4;  
MSBPR-I: threshold_i = 0.92;     
MSBPR-II: threshold_u = 0.7; threshold_c = 0.8;  
MSBPR-III: threshold_u = 0.7; threshold_c = 0.4;  	

### tmall-5k5k
fnTrainData = "./data/Tmall/tmall-5k5k-train";  
fnTestData = "./data/Tmall/tmall-5k5k-test";  
fnValidData = "./data/Tmall/tmall-5k5k-test";  
fnSimilarity_i = "./similarity_matrix/tmall/similarity_i.txt";  
fnSimilarity_u = "./similarity_matrix/tmall/similarity_u.txt";  
fnSimilarity_i_u = "./similarity_matrix/tmall/similarity_i_u.txt";  
n = 4534, m = 5000 ;  
alpha_u = alpha_v = beta_v = 0.01,  gamma = 0.01;  
threshold_p = 4;  
MSBPR-I:   threshold_i = 0.5;  
MSPBR_II:  threshold_u = 0.6, threshold_c = 0.9;  
MSBPR-III: threshold_u = 0.6, threshold_c = 0.2;  

### rsc15-5k5k
fnTrainData = "./data/rsc15/rsc15-clicks-5k5k-train";  
fnTestData = "./data/rsc15/rsc15-clicks-5k5k-test";  
fnValidData = "./data/rsc15/rsc15-clicks-5k5k-test";  
fnSimilarity_i = "./similarity_matrix/rsc15/similarity_i.txt";  
fnSimilarity_u = "./similarity_matrix/rsc15/similarity_u.txt";  
fnSimilarity_i_u = ./similarity_matrix/rsc15/similarity_i_u.txt";  
n = 4926, m = 4997 ; 
alpha_u = alpha_v = beta_v = 0.01,  gamma = 0.01;  
threshold_p = 1;  
MSBPR-I:   threshold_i = 0.98;    
MSPBR_II:  threshold_u = 0.9, threshold_c = 0.5;  
MSBPR-III: threshold_u = 0.5, threshold_c = 0.7;  

### aotm-5k5k
fnTrainData = "./data/aotm/aotm-5k5k-train";  
fnTestData = "D:/Java/MBPR/data/aotm/aotm-5k5k-test";  
fnValidData = "D:/Java/MBPR/data/aotm/aotm-5k5k-test";  
fnSimilarity_i = "D:/Java/MBPR/similarity_matrix/aotm/similarity_i.txt";  
fnSimilarity_u = "D:/Java/MBPR/similarity_matrix/aotm/similarity_u.txt";  
fnSimilarity_i_u = "D:/Java/MBPR/similarity_matrix/aotm/similarity_i_u.txt";  
public static int n = 4885 ;  
public static int m = 5000 ;  
alpha_u = alpha_v = beta_v = 0.01,  gamma = 0.01;  
threshold_p = 2;  
MSBPR-I:   threshold_i = 0.96;  
MSPBR_II:  threshold_u = 0.6, threshold_c = 0.3;  
MSBPR-III: threshold_u = 0.9, threshold_c = 0.5;  

### yelp-5k5k
fnTrainData = "./data/yelp/yelp-5k5k-train";  
fnTestData = "./data/yelp/yelp-5k5k-test";  
fnValidData = "./data/yelp/yelp-5k5k-valid";  
fnSimilarity_i = "./similarity_matrix/yelp/similarity_i.txt";  
fnSimilarity_u = "./similarity_matrix/yelp/similarity_u.txt";  
fnSimilarity_i_u = "./similarity_matrix/yelp/similarity_i_u.txt";  
public static int n = 4978 ;  
public static int m = 5000 ;  
alpha_u = alpha_v = beta_v = 0.01,  gamma = 0.01;  
threshold_p = 1;  
MSBPR-I:   threshold_i = 1;  
MSPBR_II:  threshold_u = 0.6, threshold_c = 0.5;  
MSBPR-III: threshold_u = 0.7, threshold_c = 0.7;  

### Last.fm-5k5k
fnTrainData = "./data/last.fm/last.fm-5k5k-train";  
fnTestData = "./data/last.fm/last.fm-5k5k-test";  
fnValidData = "./data/last.fm/last.fm-5k5k-valid";  
fnSimilarity_i = "./similarity_matrix/last.fm/similarity_i.txt";  
fnSimilarity_u = "./similarity_matrix/last.fm/similarity_u.txt";  
fnSimilarity_i_u = "./similarity_matrix/last.fm/similarity_i_u.txt";  
public static int n = 1640 ;   
public static int m = 12440 ;  
alpha_u = alpha_v = beta_v = 0.01,  gamma = 0.01;  
threshold_p = 5;  
MSBPR-I:   threshold_i = 0.5;  
MSPBR_II:  threshold_u = 0.9, threshold_c = 0.9;  
MSBPR-III: threshold_u = 0.5, threshold_c = 0.4;  

## selecting for similarity calculation method
(Which method to choose, mark true accordingly and the remaining two as false)  
flagsimilarity_i = true;    
flagsimilarity_u = false;    
flagsimilarity_i_u = false;   

## selecting for split method
(Select the item division method corresponding to the similarity calculation method)  
flagsplit_i = true;  
flagsplit_u = false;    
flagsplit_i_u = false;  

## Select whether the output result is test datasets or valid datasets or both output  
flagvalid = false;  
flagtest = true;  

## select the evaluation matrics, default is true  
public static boolean flagMRR = true;  
public static boolean flagMAP = true;  