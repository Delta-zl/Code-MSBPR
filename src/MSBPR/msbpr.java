package MSBPR;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;
import java.util.Map.Entry;


public class msbpr
{
	// === Configurations
	
	// ====== Example for running the whole code by setting the configurations =======
	// here we take running 30music as an example
	// the MSBPR model selects MSBPR-I;
	// the similarity calculation selects the item similarity; 
	// the item division method selects the corresponding split_i: flagsplit_i = true;
	
	// the number of latent dimensions
	public static int d = 20;
    // top k in evaluation
	public static int topK = 5;
	// regularization parameter
	public static float alpha_u = 0.01f;
    public static float alpha_v = 0.01f;
    public static float beta_v = 0.01f;
    // learning rate 
    public static float gamma = 0.01f;
    
    // number of users
    public static int n = 4586; 
    // number of items
	public static int m = 4989; 
	// number of iterations
	public static int num_iterations = 100000; 
    
    // Input threshold
    public static int threshold_c = 4;        // item popularity threshold
    public static float threshold_i = 0.92f;  // item similarity threshold
	public static float threshold_u = 0.7f;  // user similarity threshold
    public static float threshold_c2 = 0.8f; // interaction threshold

    // === Input data files
    public static String fnTrainData = "D:/Java/MBPR/data/30music/30music-5k5k-train";  // train set
    public static String fnTestData = "D:/Java/MBPR/data/30music/30music-5k5k-test";   // test set
    public static String fnValidData = "D:/Java/MBPR/data/30music/30music-5k5k-test";  // valid set
    
    // == Input similarity matrix files
    public static String fnSimilarity_i = "D:/Java/MBPR/similarity_matrix/30music/similarity_i.txt";  // item similarity
    public static String fnSimilarity_u = "D:/Java/MBPR/similarity_matrix/30music/similarity_u.txt";  // user similarity
    public static String fnSimilarity_i_u = "D:/Java/MBPR/similarity_matrix/30music/similarity_i_u.txt";  // hybrid similarity
  		
	// select the similarity calculation method(true is selected, false is not selected)
	public static boolean flagsimilarity_i = false; 
	public static boolean flagsimilarity_u = false; 
	public static boolean flagsimilarity_i_u = false;
	
	// select the split method(true is selected, false is not selected) 
	public static boolean flagsplit_i = false;
	public static boolean flagsplit_u = true;
	public static boolean flagsplit_i_u = false;
		
	// valid & test
	public static boolean flagvalid = false;
	public static boolean flagtest = true;
	
	// === Evaluation
	public static boolean flagMRR = true;
	public static boolean flagMAP = true;
	
	// average value
	public static float Pre_ave = 0;
	public static float Rec_ave = 0;
	public static float F1_ave = 0;
	public static float NDCG_ave = 0;
	public static float MRR_ave = 0;
	public static float MAP_ave = 0;

	// === Data
	public static HashMap<Integer, HashMap<Integer, Integer>> Data = new HashMap<Integer, HashMap<Integer,Integer>>();
	public static HashMap<Integer, HashMap<Integer, Integer>> DataItem2User = new HashMap<Integer, HashMap<Integer,Integer>>(); 

    // === training data
    public static HashMap<Integer, HashMap<Integer, Integer>> TrainData = new HashMap<Integer, HashMap<Integer,Integer>>(); 
    public static HashMap<Integer, HashMap<Integer, Integer>> TrainDataItem2User = new HashMap<Integer, HashMap<Integer,Integer>>(); 

	// === Iu_p data, potentially liked items
	public static HashMap<Integer, List<Integer>> Iu_p_data = new HashMap<Integer, List<Integer>>();
	// === Iu_c data, uncertainly liked items
	public static HashMap<Integer, List<Integer>> Iu_c_data = new HashMap<Integer, List<Integer>>();
	// === Iu_j data, disliked items
	public static HashMap<Integer, List<Integer>> Iu_j_data = new HashMap<Integer, List<Integer>>();
	
    // === test data
    public static HashMap<Integer, HashSet<Integer>> TestData = new HashMap<Integer, HashSet<Integer>>(); 
    // === validation data
    public static HashMap<Integer, HashSet<Integer>> ValidData = new HashMap<Integer, HashSet<Integer>>();
        
    // === whole item set
    public static HashSet<Integer> ItemSetWhole = new HashSet<Integer>();  
    // === whole user set
    public static HashSet<Integer> UserSetWhole = new HashSet<Integer>();
    
    // === The item that most users may dislike
    public static HashSet<Integer> ItemSetDislike = new HashSet<Integer>();
    
    // === some statistics, start from index "0"
    public static int[] itemRatingNumTrain; 
    
    // === model parameters to learn, start from index "0"
    public static float[][] U;
    public static float[][] V;
    public static float[] biasV;  // bias of item
        
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    public static void main(String[] args) throws Exception
    {	
    	// === Print the configurations
		System.out.println("-------" + "MSBPR" + "-------");
		
		System.out.println("（1）Parameter settings: ");
    	System.out.println("d: " + Integer.toString(d));
    	System.out.println("alpha_u: " + Float.toString(alpha_u));
    	System.out.println("alpha_v: " + Float.toString(alpha_v));
    	System.out.println("beta_v: " + Float.toString(beta_v));
    	System.out.println("gamma: " + Float.toString(gamma));
    	System.out.println("threshold_c: " + Integer.toString(threshold_c));
    	System.out.println("threshold_i: " + Float.toString(threshold_i));
		System.out.println("threshold_u: " + Float.toString(threshold_u));
    	System.out.println("threshold_c2: " + Float.toString(threshold_c2));
    	System.out.println("n: " + Integer.toString(n));
    	System.out.println("m: " + Integer.toString(m));
    	System.out.println("num_iterations: " + Integer.toString(num_iterations));
    	System.out.println("topK: " + Integer.toString(topK));
    	
    	System.out.println("（2）Selection of dataset: ");
    	System.out.println("fnTestData: " + fnTestData);
    	System.out.println("fnValidData: " + fnValidData);
    	System.out.println("fnValidData: " + fnValidData);

    	System.out.println("（3）Selection of similarity calculation: ");   	
    	System.out.println("flagsimilarity_i: " + Boolean.toString(flagsimilarity_i));
    	System.out.println("flagsimilarity_u: " + Boolean.toString(flagsimilarity_u));
    	System.out.println("flagsimilarity_i_u: " + Boolean.toString(flagsimilarity_i_u));
    	System.out.println("fnSimilarity_i: " + fnSimilarity_i);
    	System.out.println("fnSimilarity_u: " + fnSimilarity_u);
    	System.out.println("fnSimilarity_i_u: " + fnSimilarity_i_u);
    	
    	System.out.println("（4）Selection of division method: "); 
    	System.out.println("flagsplit_i: " + Boolean.toString(flagsplit_i));
    	System.out.println("flagsplit_u: " + Boolean.toString(flagsplit_u));
		System.out.println("flagsplit_i_u: " + Boolean.toString(flagsplit_i_u));
	
    	
    	System.out.println("（5）Selection of evaluation indicators: ");
    	System.out.println("flagMRR: " + Boolean.toString(flagMRR));
    	System.out.println("flagMAP: " + Boolean.toString(flagMAP));    	
    	// --------------------------------------------------------------
        itemRatingNumTrain = new int[m + 1]; // start from index "1"

		// === Locate memory for the data structure of the model parameters
        U = new float[n+1][d];
        V = new float[m+1][d];
        biasV = new float[m+1];  // bias of item
    	// --------------------------------------------------------------     	
		
				
        // === Step 1: Read data
    	long TIME_START_READ_DATA = System.currentTimeMillis();
    	readDataTrainTestValid();
    	long TIME_FINISH_READ_DATA = System.currentTimeMillis();
    	System.out.println("Time (read data):" + 
    				Float.toString((TIME_FINISH_READ_DATA-TIME_START_READ_DATA)/1000F)
    				+ "s");
    	// ------------------------------   	
  
    	// ------------------------------
    	// === Step 2: Initialization of U, V, bias
    	long TIME_START_INITIALIZATION = System.currentTimeMillis();
    	initialize();
    	long TIME_FINISH_INITIALIZATION = System.currentTimeMillis();
    	System.out.println("Time (initialization):" + 
    				Float.toString((TIME_FINISH_INITIALIZATION-TIME_START_INITIALIZATION)/1000F)
    				+ "s");
    	// ------------------------------

        // === Step 3: Caculation of Similarity
        long TIME_START_SIM = System.currentTimeMillis();
        similarity_i();
        similarity_u();
   	    similarity_i_u();
        long TIME_FINISH_SIM = System.currentTimeMillis();
        System.out.println("Time (similarity):" +
               Float.toString((TIME_FINISH_SIM-TIME_START_SIM)/1000F)
               + "s");
        // ------------------------------
    	
		// ------------------------------
		// === Step 4: Spliting	items	
		long TIME_START_SPLIT = System.currentTimeMillis();
		split_i();
		split_u();
		split_i_u();
		long TIME_FINISH_SPLIT = System.currentTimeMillis();
		System.out.println("Time (spliting):" +
				Float.toString((TIME_FINISH_SPLIT-TIME_START_SPLIT)/1000F)
				+ "s");
        // ------------------------------
		
		// ------------------------------
		// === Step 5: Training
		long TIME_START_TRAIN = System.currentTimeMillis();
		train();
		long TIME_FINISH_TRAIN = System.currentTimeMillis();
		System.out.println("Time (training):" +
				Float.toString((TIME_FINISH_TRAIN - TIME_START_TRAIN) / 1000F) + "s");
    	// ------------------------------

    	// ------------------------------
    	// === Step 6: Prediction and Evaluation
    	if (flagvalid)
    	{	
    		System.out.println("--- valid ---");
    		long TIME_START_VALID = System.currentTimeMillis();
	    	testRanking(ValidData);
	    	long TIME_FINISH_VALID = System.currentTimeMillis();
	    	System.out.println("Time (validation):" +
	    				Float.toString((TIME_FINISH_VALID-TIME_START_VALID)/1000F)
	    				+ "s");
    	}
		// === Step 7: TestRanking
    	if (flagtest)
    	{
    		System.out.println("--- test ---");
    		long TIME_START_TEST = System.currentTimeMillis();
    		testRanking(TestData);
    		long TIME_FINISH_TEST = System.currentTimeMillis();
    		System.out.println("Time (test):" +
    					Float.toString((TIME_FINISH_TEST-TIME_START_TEST)/1000F)
    					+ "s");
    	}
    }
          
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    public static void readDataTrainTestValid() throws Exception
    {	
    	// ---------------------TrainData------------------------
    	BufferedReader br = new BufferedReader(new FileReader(fnTrainData));
    	String line = null;    			
    	while ((line = br.readLine())!=null)
    	{
    		String[] terms = line.split("\t");
    		int userID = Integer.parseInt(terms[0]);    		
    		int itemID = Integer.parseInt(terms[1]);
    		int count = Integer.parseInt(terms[2]);
    		
    		// --- add to the whole item set
    		ItemSetWhole.add(itemID);
		    			
    		// --- add to the whole user set
    		UserSetWhole.add(userID);
    		
		    // --- statistics, used to calculate the interacitions of each item for all users
		    itemRatingNumTrain[itemID] += count;  
    		   		
			// TrainData: user->items   
			if(TrainData.containsKey(userID))
			{
				HashMap<Integer, Integer> itemset = TrainData.get(userID);    	
	    		itemset.put(itemID, count);
	    		TrainData.put(userID, itemset);
	    	}
	    	else
	    	{
				HashMap<Integer, Integer> itemset = new HashMap<Integer, Integer>();
				itemset.put(itemID, count);
				TrainData.put(userID, itemset);
	    	}

			// Data: user->items   
			if(Data.containsKey(userID))
			{
				HashMap<Integer, Integer> itemset = Data.get(userID);
				itemset.put(itemID, count);
				Data.put(userID, itemset);
			}
			else
			{
				HashMap<Integer, Integer> itemset = new HashMap<Integer, Integer>();
				itemset.put(itemID, count);
				Data.put(userID, itemset);
			}

			// TrainDataItem2User: item->users   
			if(TrainDataItem2User.containsKey(itemID))
    	    {
    	    	HashMap<Integer, Integer> userSet = TrainDataItem2User.get(itemID);
    	    	if (userSet.size()<10000)
    	    	{
        	   		userSet.put(userID, count);
        	   		TrainDataItem2User.put(itemID, userSet);
    	    	}
    	    }
    	    else
    	    {
    	    	HashMap<Integer, Integer> userSet = new HashMap<Integer, Integer>();
    	    	userSet.put(userID, count);
    	    	TrainDataItem2User.put(itemID, userSet);
    	    }
			// DataItem2User: item->users   
			if(DataItem2User.containsKey(itemID))
    	    {
    	    	HashMap<Integer, Integer> userSet = DataItem2User.get(itemID);
    	    	if (userSet.size()<10000)
    	    	{
        	   		userSet.put(userID, count);
        	   		DataItem2User.put(itemID, userSet);
    	    	}
    	    }
    	    else
    	    {
    	    	HashMap<Integer, Integer> userSet = new HashMap<Integer, Integer>();
    	    	userSet.put(userID, count);
    	    	DataItem2User.put(itemID, userSet);
    	    }
    	}
    	br.close();
    	// ======================================================
		// ==================== Itemdislike creation===================
		for(Integer i:TrainDataItem2User.keySet())
		{
			if(TrainDataItem2User.get(i).size() <= threshold_c)
			{
				ItemSetDislike.add(i);
			}
		}
		System.out.println("size_itemdislike: " + ItemSetDislike.size());
		// ------------------------------------------------------

    	// ---------------------TestData-------------------------
    	if (fnTestData.length()>0) {
			br = new BufferedReader(new FileReader(fnTestData));
			String line1 = null;
			while ((line1 = br.readLine()) != null) {
				String[] terms = line1.split("\t");
				int userID = Integer.parseInt(terms[0]);
				int itemID = Integer.parseInt(terms[1]);
				int count = Integer.parseInt(terms[2]);

				// --- add to the whole item set
				ItemSetWhole.add(itemID);

				// --- add to the whole user set
				UserSetWhole.add(userID);

				// --- test data
				if (TestData.containsKey(userID)) {
					HashSet<Integer> itemSet = TestData.get(userID);
					itemSet.add(itemID);
					TestData.put(userID, itemSet);
				} else {
					HashSet<Integer> itemSet = new HashSet<Integer>();
					itemSet.add(itemID);
					TestData.put(userID, itemSet);
				}

				// Data: user->items   
				if(Data.containsKey(userID))
				{
					HashMap<Integer, Integer> itemset = Data.get(userID);
					itemset.put(itemID, count);
					Data.put(userID, itemset);
				}
				else
				{
					HashMap<Integer, Integer> itemset = new HashMap<Integer, Integer>();
					itemset.put(itemID, count);
					Data.put(userID, itemset);
				}
				// DataItem2User: item->users   
				if(DataItem2User.containsKey(itemID))
	    	    {
	    	    	HashMap<Integer, Integer> userSet = DataItem2User.get(itemID);
	    	    	if (userSet.size()<10000)
	    	    	{
	        	   		userSet.put(userID, count);
	        	   		DataItem2User.put(itemID, userSet);
	    	    	}
	    	    }
	    	    else
	    	    {
	    	    	HashMap<Integer, Integer> userSet = new HashMap<Integer, Integer>();
	    	    	userSet.put(userID, count);
	    	    	DataItem2User.put(itemID, userSet);
	    	    }
			}
			br.close();
		}
		// =======================================================
   	
    	// ----------------------ValidData------------------------
    	if (fnValidData.length()>0)
    	{
	    	br = new BufferedReader(new FileReader(fnValidData));
	    	line = null;
	    	while ((line = br.readLine())!=null)
	    	{
	    		String[] terms = line.split("\t");
	    		int userID = Integer.parseInt(terms[0]);
	    		int itemID = Integer.parseInt(terms[1]);  
	    		int count = Integer.parseInt(terms[2]);
	    	
	    		// --- add to the whole item set
				ItemSetWhole.add(itemID);

				// --- add to the whole user set
				UserSetWhole.add(userID);
				
				// --- validation data
				if(ValidData.containsKey(userID))
	    		{
	    			HashSet<Integer> itemSet = ValidData.get(userID);
	    			itemSet.add(itemID);
	    			ValidData.put(userID, itemSet);
	    		}
	    		else
	    		{
	    			HashSet<Integer> itemSet = new HashSet<Integer>();
	    			itemSet.add(itemID);
	    			ValidData.put(userID, itemSet);
	    		}
				// Data: user->items   
				if(Data.containsKey(userID))
				{
					HashMap<Integer, Integer> itemset = Data.get(userID);
					itemset.put(itemID, count);
					Data.put(userID, itemset);
				}
				else
				{
					HashMap<Integer, Integer> itemset = new HashMap<Integer, Integer>();
					itemset.put(itemID, count);
					Data.put(userID, itemset);
				}
				// DataItem2User: item->users   
				if(DataItem2User.containsKey(itemID))
	    	    {
	    	    	HashMap<Integer, Integer> userSet = DataItem2User.get(itemID);
	    	    	if (userSet.size()<10000)
	    	    	{
	        	   		userSet.put(userID, count);
	        	   		DataItem2User.put(itemID, userSet);
	    	    	}
	    	    }
	    	    else
	    	    {
	    	    	HashMap<Integer, Integer> userSet = new HashMap<Integer, Integer>();
	    	    	userSet.put(userID, count);
	    	    	DataItem2User.put(itemID, userSet);
	    	    }

	    	}
	    	br.close();
    	}
    	// ----------------------------------------------------
    	
    }
    

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    public static void initialize()
    {
    	// ======================================================
    	Random r = new Random(1);
    	// --- initialization of U and V
    	for (int u=1; u<n+1; u++)
    	{
    		for (int f=0; f<d; f++)
    		{
    			U[u][f] = (float) ( Math.sqrt(0.01) * r.nextGaussian());
    		}
    	}
    	//
    	for (int i=1; i<m+1; i++)
    	{
    		for (int f=0; f<d; f++)
    		{
    			V[i][f] = (float) ( Math.sqrt(0.01) * r.nextGaussian());
    		}
    	}
    	// ======================================================

    	// ======================================================
    	// --- initialization of \mu and b_i
    	float g_avg = 0;
    	for (int i=1; i<m+1; i++)
    	{
    		g_avg += itemRatingNumTrain[i];
    	}
    	g_avg = g_avg/n/m;  
    	System.out.println( "The global average rating:" + Float.toString(g_avg) );

    	// --- biasV[i] represents the popularity of the item i, which is initialized to [0,1]
    	for (int i=1; i<m+1; i++)
    	{
    		 biasV[i]= (float) itemRatingNumTrain[i] / n - g_avg;
    	}
    	// $ \mu = \sum_{u,i} y_{ui} /n/m $
    	// $ b_i = \sum_{u=1}^n (y_{ui} - \mu) / n $
    	// ======================================================
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    public static void similarity_i() {
    	/**
    	 * The similarity matrix between items is calculated in advance and stored in the file, 
    	 * which saves time and facilitates later calling and parameter adjustment.
    	 * */
    	if(flagsimilarity_i) 
    	{	
	    	System.out.println("similarity:  similarity_item");  
	    	double [][] S = new double [m+1][m+1];
	     	double fenzi = 0;
	    	double fenmu1 = 0;
	    	double fenmu2 = 0;
	    	List<Integer> list1 = new ArrayList<Integer>(); 
	    	List<Integer> list2 = new ArrayList<Integer>(); 
	    	
	    	// -----------------------------------
	    	for(int i1 = 1; i1 <= m; i1++)
	    	{
				if (!DataItem2User.containsKey(i1)){    
					continue;
				}
	    		List<Integer> userset_i1 = new ArrayList<Integer>(DataItem2User.get(i1).keySet());
	    		
	    		for(int i2 = i1+1; i2 <= m; i2++) 
	    		{
    				list1.clear();
    				list2.clear();    				
	    			fenzi = fenmu1 = fenmu2 = 0;
					if (!DataItem2User.containsKey(i2)){   
						continue;
					}
        			List<Integer> userset_i2 = new ArrayList<Integer>(DataItem2User.get(i2).keySet());

        			for(Integer u:userset_i1) {
        				if(userset_i2.contains(u)) {

        					fenzi += Data.get(u).get(i1) * Data.get(u).get(i2);
        					fenmu1 += Data.get(u).get(i1) * Data.get(u).get(i1);
        					fenmu2 += Data.get(u).get(i2) * Data.get(u).get(i2);
        				}
        			}
        			if (fenzi == 0)  
        			{
        				S[i1][i2] = 0;
        			}
        			else
        			{
        				S[i1][i2] = fenzi / (Math.sqrt(fenmu1) * Math.sqrt(fenmu2));
        			}
        	
	   			System.out.println("S(i)" + i1 + "-" + i2 + ":  " + S[i1][i2]);
	    		}
    	}
    	
    	try {    		
			System.out.println("writing(similarity_i)..........");
			
			PrintStream out = System.out;

			PrintStream ps = new PrintStream(fnSimilarity_i);
			System.setOut(ps);
			
			for(int i=0; i<S.length; i++) 
			{
				for(int j=0; j<S[0].length; j++) 
				{
					if(j==i) {
						S[i][j]=1;
						System.out.print(S[i][j] + "\t");
					}
					else if(j>i){
						System.out.print(S[i][j] + "\t");
					}
					else if(j<i) {
						S[i][j] = S[j][i];
						System.out.print(S[i][j] + "\t");
					}
				}
				System.out.println();
			}
			
			System.setOut(out);
			System.out.println("file (similarity_i) record complete, please check the corresponding file.");
    		
    	}catch(FileNotFoundException e) {
    		e.printStackTrace();
    	}
      }
    }
   
    // =========================================================
    public static void split_i() throws Exception
    {
    	if(flagsplit_i) 
    	{

    	System.out.println("split:  split_i");
    	double [][] S = new double [m+1][m+1];
    	int n_ = 0;

    	BufferedReader br = new BufferedReader(new FileReader(fnSimilarity_i));
    	String line = null;    			
    	while((line = br.readLine())!= null) {
    		
    		String[] terms = line.split("\t");
    		
    		for(int i=0; i<terms.length; i++) {
    			
    			S[n_][i] = Double.parseDouble(terms[i]);
    		}
    		
    		n_++;
    	}   	  	
    	System.out.println("Similarity matrix size: " + "row: " + S.length + " col: " + S[0].length);	
    	
    	// -------------------- split item --------------------------
        for (int u = 1; u <= n; u++)
        {
            // ------------------------Iu_i---------------------------
            List<Integer> Iu_i = new ArrayList<Integer>(TrainData.get(u).keySet());    // === Iu_i
            int Size_Iu_i = Iu_i.size();                                               
            List<Integer> Iu_i_count = new ArrayList<Integer>(TrainData.get(u).values()); 
            int sum = 0;
            for(Integer c:Iu_i_count) {
                sum += c;                                                           
            }

            // ------------------------Iu_j---------------------------
            List<Integer> Iu_j = new ArrayList<Integer>(ItemSetDislike);               
            Iu_j.removeAll(Iu_i);                                                      // === Iu_j
            Iu_j_data.put(u, Iu_j);
            int Size_Iu_j = Iu_j.size();

            // ----------------------Iu_*------------------
            List<Integer> Iu_unknow = new ArrayList<Integer>(ItemSetWhole);            
            Iu_unknow.removeAll(Iu_i);
            Iu_unknow.removeAll(Iu_j);                                                 // === Iu_*
            int Size_Iu_unknow = Iu_unknow.size();                                    

            // -------------------------similarity--------------------
            // --- caculate the similarity_i with weight
			HashMap<Integer, Double> similarity = new HashMap<Integer, Double>();  
			similarity.clear();
			
            for(Integer i1:Iu_i)
            {
                double weight = TrainData.get(u).get(i1);
                weight = weight/sum;  

                for(Integer i2:Iu_unknow)
                {
                	double cosine = S[i1][i2] * weight;

                    if(similarity.containsKey(i2))
                    {
                        Double cos = similarity.get(i2);
                        cos += cosine;
                        similarity.put(i2, cos);
                    }
                    else
                    {
                        similarity.put(i2, cosine);
                    }

                }
            }

            // -----------------------sort----------------------
            List<Map.Entry<Integer,Double>> sorted_similarity =
                    new ArrayList<Map.Entry<Integer,Double>>(similarity.entrySet());
            Collections.sort(sorted_similarity, new Comparator<Map.Entry<Integer,Double>>()
            {
                public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
                {
                    return o2.getValue().compareTo( o1.getValue() );
                }
            });	    	

            // -----------------------Iu_p---------------------------
            List<Integer> Iu_p = new ArrayList<Integer>();
            for(Entry<Integer, Double> i4:sorted_similarity) {
                if(i4.getValue() >= threshold_i)
                {
                    Iu_p.add(i4.getKey());                                   // === Iu_p
                }
            }
            Iu_p_data.put(u, Iu_p);

            // -----------------------Iu_c--------------------------
            List<Integer> Iu_c = new ArrayList<Integer>(Iu_unknow);
            Iu_c.removeAll(Iu_p);                                           // === Iu_c
            Iu_c_data.put(u, Iu_c);
            
            // ----------
//            System.out.println("Iu_i:" + Size_Iu_i + Iu_i);
//            System.out.println("Iu_p: " + Iu_p.size() + Iu_p);
//            System.out.println("Iu_c: " + Iu_c.size());
//            System.out.println("Iu_j: " + Size_Iu_j + Iu_j);

        }
      }
    }
     	
	// ========================================================
    public static void similarity_u() {
    	
    	/**
    	 * The similarity matrix between items is calculated in advance and stored in the file, 
    	 * which saves time and facilitates later calling and parameter adjustment.
    	 * */
    	if(flagsimilarity_u)
    	{
    	System.out.println("similarity:  similarity_user");
    	double [][] S = new double [n+1][n+1];
    	double fenzi = 0;
    	double fenmu1 = 0;
    	double fenmu2 = 0;
    	List<Integer> list1 = new ArrayList<Integer>(); 
        List<Integer> list2 = new ArrayList<Integer>(); 

    	for (int u1 = 1; u1 <= n; u1++)
    	{
			if (!Data.containsKey(u1)){    
				continue;
			}
    		List<Integer> Iu1_i = new ArrayList<Integer>(Data.get(u1).keySet());     

    		for (int u2 = u1+1; u2 <= n; u2++)
    		{   		
    			list1.clear();
    			list2.clear();   		
    			fenzi = fenmu1 = fenmu2 = 0;
				if (!Data.containsKey(u2)){   
					continue;
				}
    			List<Integer> Iu2_i = new ArrayList<Integer>(Data.get(u2).keySet());   
    			for(Integer i:Iu1_i)
    			{
    				if(Iu2_i.contains(i))
    				{                                                  
    					fenzi += Data.get(u1).get(i) * Data.get(u2).get(i);
    					fenmu1 += Data.get(u1).get(i) * Data.get(u1).get(i);
    					fenmu2 += Data.get(u2).get(i) * Data.get(u2).get(i);
    				}
    			}
    			if (fenzi == 0)  
    			{
    				S[u1][u2] = 0;
    			}
    			else            
    			{
    				S[u1][u2] = fenzi / (Math.sqrt(fenmu1) * Math.sqrt(fenmu2));
    			}

    			System.out.println("S(u)" + u1 + "-" + u2 + ": " + S[u1][u2]);

    		}
    	}
    	
    	// -------------- Write the similarity matrix to a file ---------------
    	try {
    		
			System.out.println("writing(similarity_u)..........");
			
			PrintStream out = System.out;

			PrintStream ps = new PrintStream(fnSimilarity_u);
			System.setOut(ps);
		
			for(int i=0; i<S.length; i++) 
			{
				for(int j=0; j<S[0].length; j++) 
				{
					if(j==i) {
						S[i][j]=0;
						System.out.print(S[i][j] + "\t");
					}
					else if(j>i){
						System.out.print(S[i][j] + "\t");
					}
					else if(j<i) {
						S[i][j] = S[j][i];
						System.out.print(S[i][j] + "\t");
					}
				}
				System.out.println();
			}
			
			System.setOut(out);
			System.out.println("file (similarity_u) record complete, please check the corresponding file.");
    		
    	}catch(FileNotFoundException e) {
    		e.printStackTrace();
    	}
      }
    }
    // ========================================================

	public static void split_u () throws Exception
    {
    	if(flagsplit_u) 
    	{
    	System.out.println("split:  split_u");

    	double [][] S = new double [n+1][n+1];
    	int n_ = 0;
    	
		BufferedReader br = new BufferedReader(new FileReader(fnSimilarity_u));
    	String line = null;    			
    	while((line = br.readLine())!= null) {
    		
    		String[] terms = line.split("\t");
    		
    		for(int i=0; i<terms.length; i++) {
    			
    			S[n_][i] = Double.parseDouble(terms[i]);
    		}
    		
    		n_++;
    	}br.close();   	    	
    	System.out.println("User similarity matrix size：" + "row: " + S.length + " col: " + S[0].length);	
			   	
    	// -------------- split items ----------
		HashMap<Integer, Double> u_s = new HashMap<Integer, Double>();  
		List<Integer> similar_user = new ArrayList<Integer>();   

        for (int u = 1; u <= n; u++)
        {
			similar_user.clear();			
			u_s.clear();

			for(int u1 = 1; u1 <= n; u1++)
			{
				u_s.put(u1, S[u][u1]);   
			}

			// ------------- Sort by user similarity and filter out similar users ----------
			List<Map.Entry<Integer,Double>> sorted_user = new ArrayList<Map.Entry<Integer,Double>>(u_s.entrySet());
			Collections.sort(sorted_user, new Comparator<Map.Entry<Integer,Double>>()
			{
				public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
				{
					return o2.getValue().compareTo( o1.getValue() ); 
				}
			});

			int sum_num = 0;
            for(Entry<Integer, Double> i4:sorted_user) {
                if(i4.getValue() >= threshold_u)
                {
					sum_num += 1;
                    similar_user.add(i4.getKey());                                   
                }
            }
			// ------------------------------------------------------

            HashMap<Integer, Double> r_unknow_item = new HashMap<Integer, Double>();     

            // ------------------------Iu_i---------------------------
            List<Integer> Iu_i = new ArrayList<Integer>(TrainData.get(u).keySet());    // === Iu_i
            int Size_Iu_i = Iu_i.size();                                               
            List<Integer> Iu_i_count = new ArrayList<Integer>(TrainData.get(u).values()); 
            int sum = 0;
            for(Integer c:Iu_i_count) {
                sum += c;                                                            
            }
            
            // ------------------------Iu_j---------------------------
            List<Integer> Iu_j = new ArrayList<Integer>(ItemSetDislike);               
            Iu_j.removeAll(Iu_i);                                                      // === Iu_j
            Iu_j_data.put(u, Iu_j);
            int Size_Iu_j = Iu_j.size();

            // ----------------------Iu_*------------------
            List<Integer> Iu_unknow = new ArrayList<Integer>(ItemSetWhole);           
            Iu_unknow.removeAll(Iu_i);
            Iu_unknow.removeAll(Iu_j);                                                 // === Iu_*
            int Size_Iu_unknow = Iu_unknow.size();                     

            // -------------------------Predict the number of interactions--------------------
			for(Integer i:Iu_unknow)
			{
				double r=0, s_sum=0;
				for(Integer uu:similar_user)      
				{
					s_sum += S[u][uu];            
					if(TrainData.get(uu).get(i) == null)
					{
						r += 0;  
					}
					else
					{
						r += S[u][uu] * TrainData.get(uu).get(i);  
					}
				}
				r = r/s_sum;
            	r_unknow_item.put(i, r);
			}

            // -----------------------sort----------------------
            List<Map.Entry<Integer,Double>> sorted_r_unknow_item =
                    new ArrayList<Map.Entry<Integer,Double>>(r_unknow_item.entrySet());
            Collections.sort(sorted_r_unknow_item, new Comparator<Map.Entry<Integer,Double>>()
            {
                public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
                {
                    return o2.getValue().compareTo( o1.getValue() ); 
                }
            });

            // -----------------------Iu_p---------------------------
            List<Integer> Iu_p = new ArrayList<Integer>();
            for(Entry<Integer, Double> i4:sorted_r_unknow_item) {
                if(i4.getValue() >= threshold_c2)
                {
                    Iu_p.add(i4.getKey());                                   // === Iu_p
                }
            }
            Iu_p_data.put(u, Iu_p);

            // -----------------------Iu_c--------------------------
            List<Integer> Iu_c = new ArrayList<Integer>(Iu_unknow);
            Iu_c.removeAll(Iu_p);                                           // === Iu_c
            Iu_c_data.put(u, Iu_c);
        }
      }
    }

	// ========================================================	
    public static void similarity_i_u() throws Exception{

        /* *
         ** The similarity matrix between items is calculated in advance and stored in the file, 
    	 * which saves time and facilitates later calling and parameter adjustment.
         * */
    	if(flagsimilarity_i_u) 
    	{

        double [][] S_i = new double [m+1][m+1]; 
        double [][] S_u = new double [n+1][n+1]; 
        double a, b ;  
        double fenzi = 0;
        int n_ = 0;

        // ------- Read out the item similarity matrix data which use it later ----------
        BufferedReader br = new BufferedReader(new FileReader(fnSimilarity_i));
        String line = null;
        while((line = br.readLine())!= null) {

            String[] terms = line.split("\t");

            for(int i=0; i<terms.length; i++) {

                S_i[n_][i] = Double.parseDouble(terms[i]);
            }

            n_++;
        }

        System.out.println("size of item similarity matrix：" + "row: " + S_i.length + " col: " + S_i[0].length);

        // ------------------------- caculate the hybrid similarity -------------------------
        for(int u1 = 1; u1 <= n; u1++){
			if (!Data.containsKey(u1)){
				continue;
			}
            List<Integer> Iu1_items = new ArrayList<Integer>(Data.get(u1).keySet());
            List<Integer> Iu1_items_values = new ArrayList<Integer>(Data.get(u1).values());
            int num_u1_items = Iu1_items.size();
            a = sum(Iu1_items_values);  
            
            for(int u2 = u1 + 1; u2 <= n; u2++){        	
            	fenzi = 0;
				if (!Data.containsKey(u2)){    
					continue;
				}
                List<Integer> Iu2_items = new ArrayList<Integer>(Data.get(u2).keySet());
                List<Integer> Iu2_items_values = new ArrayList<Integer>(Data.get(u2).values());
                int num_u2_items = Iu2_items.size();
                b = sum(Iu2_items_values);  
                S_u[u1][u2] = (Math.sqrt(a*num_u2_items)) * (Math.sqrt(b*num_u1_items));  

                for (Integer i1:Iu1_items){
                	                
                    for (Integer i2:Iu2_items){

                        fenzi += Data.get(u1).get(i1) * Data.get(u2).get(i2) * S_i[i1][i2];                        
                    }
                }
                S_u[u1][u2] = fenzi / S_u[u1][u2];                                
                System.out.println("S(i_u)" + u1 + "-" + u2 + ": " + S_u[u1][u2]);
                
            }
        }
    	try {
    		
			System.out.println("writing(similarity_i_u)..........");
			
			PrintStream out = System.out;
			PrintStream ps = new PrintStream(fnSimilarity_i_u);
			System.setOut(ps);
			
			for(int i=0; i<S_u.length; i++) 
			{
				for(int j=0; j<S_u[0].length; j++) 
				{
					if(j==i) {
						S_u[i][j]=0;
						System.out.print(S_u[i][j] + "\t");
					}
					else if(j>i){
						System.out.print(S_u[i][j] + "\t");
					}
					else if(j<i) {
						S_u[i][j] = S_u[j][i];
						System.out.print(S_u[i][j] + "\t");
					}
				}
				System.out.println();
			}		
			System.setOut(out);
			System.out.println("file(similarity_i_u)record complete, please check the corresponding file.");
    		
    	}catch(FileNotFoundException e) {
    		e.printStackTrace();
    	}
      }
    }
    
	// ========================================================
	public static void split_i_u () throws Exception
    {
    	if(flagsplit_i_u) 
    	{
    	System.out.println("split:  split_i_u");

    	double [][] S = new double [n+1][n+1];
    	int n_ = 0;
    	
		BufferedReader br = new BufferedReader(new FileReader(fnSimilarity_i_u));
    	String line = null;    			
    	while((line = br.readLine())!= null) {
    		
    		String[] terms = line.split("\t");
    		
    		for(int i=0; i<terms.length; i++) {
    			
    			S[n_][i] = Double.parseDouble(terms[i]);
    		}
    		
    		n_++;
    	}br.close();   	
			   	
    	// -------------- Predict the number of interactions based on similarity and divide the dataset ----------
		HashMap<Integer, Double> u_s = new HashMap<Integer, Double>(); 
		List<Integer> similar_user = new ArrayList<Integer>();   

        for (int u = 1; u <= n; u++)
        {
			similar_user.clear();			
			u_s.clear();

			for(int u1 = 1; u1 <= n; u1++)
			{
				u_s.put(u1, S[u][u1]);   
			}

			// ------------- Sort by user similarity and filter out similar users --------
			List<Map.Entry<Integer,Double>> sorted_user = new ArrayList<Map.Entry<Integer,Double>>(u_s.entrySet());
			Collections.sort(sorted_user, new Comparator<Map.Entry<Integer,Double>>()
			{
				public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
				{
					return o2.getValue().compareTo( o1.getValue() ); 
				}
			});

			int sum_num = 0;
            for(Entry<Integer, Double> i4:sorted_user) {
                if(i4.getValue() >= threshold_u)
                {
					sum_num += 1;
                    similar_user.add(i4.getKey());                                   
                }
            }
			// ------------------------------------------------------

            HashMap<Integer, Double> r_unknow_item = new HashMap<Integer, Double>();     

            // ------------------------Iu_i---------------------------
            List<Integer> Iu_i = new ArrayList<Integer>(TrainData.get(u).keySet());    // === Iu_i
            int Size_Iu_i = Iu_i.size();                                              
            List<Integer> Iu_i_count = new ArrayList<Integer>(TrainData.get(u).values()); 
            int sum = 0;
            for(Integer c:Iu_i_count) {
                sum += c;                                                          
            }

            // ------------------------Iu_j---------------------------
            List<Integer> Iu_j = new ArrayList<Integer>(ItemSetDislike);               
            Iu_j.removeAll(Iu_i);                                                      // === Iu_j
            Iu_j_data.put(u, Iu_j);
            int Size_Iu_j = Iu_j.size();

            // ----------------------Iu_*(Iu_unknow)------------------
            List<Integer> Iu_unknow = new ArrayList<Integer>(ItemSetWhole);           
            Iu_unknow.removeAll(Iu_i);
            Iu_unknow.removeAll(Iu_j);                                                 // === Iu_*
            int Size_Iu_unknow = Iu_unknow.size();                                    

            // -------------------------Predict the number of interactions--------------------
			for(Integer i:Iu_unknow)
			{
				double r=0, s_sum=0;
				for(Integer uu:similar_user)      
				{
					s_sum += S[u][uu];            
					if(Data.get(uu).get(i) == null)
					{
						r += 0;  
					}
					else
					{
						r += S[u][uu] * Data.get(uu).get(i);  
					}
				}
				r = r/s_sum;
            	r_unknow_item.put(i, r);
			}

            // -----------------------sort----------------------
            List<Map.Entry<Integer,Double>> sorted_r_unknow_item =
                    new ArrayList<Map.Entry<Integer,Double>>(r_unknow_item.entrySet());
            Collections.sort(sorted_r_unknow_item, new Comparator<Map.Entry<Integer,Double>>()
            {
                public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
                {
                    return o2.getValue().compareTo( o1.getValue() ); 
                }
            });

            // -----------------------Iu_p---------------------------
            List<Integer> Iu_p = new ArrayList<Integer>();
            for(Entry<Integer, Double> i4:sorted_r_unknow_item) {
                if(i4.getValue() >= threshold_c2)
                {
                    Iu_p.add(i4.getKey());                                   // === Iu_p
                }
            }
            Iu_p_data.put(u, Iu_p);

            // -----------------------Iu_c--------------------------
            List<Integer> Iu_c = new ArrayList<Integer>(Iu_unknow);
            Iu_c.removeAll(Iu_p);                                           // === Iu_c
            Iu_c_data.put(u, Iu_c);

        }
      }
    }

	//  ========================================================
	public static void train()
	{
		float pre_total_loss = 0f;
		float total_loss = 0f;
        for (int iter = 1; iter <= num_iterations; iter++)
        {
		for (int iter_rand = 1; iter_rand <= n; iter_rand++)
		{
			// ------------------------------- u --------------------------------
			// --- randomly sample a user $u$, Math.random(): [0.0, 1.0)
			int u = (int) Math.floor(Math.random() * n) + 1; 
			if (!TrainData.containsKey(u))
				continue;
			// ==================================================================
			// ------------------------------ Iu_i -------------------------------
			List<Integer> Iu_i = new ArrayList<Integer>(TrainData.get(u).keySet());    // === Iu_i
			int Size_Iu_i = Iu_i.size();
			// --- randomly sample an item $i$, Math.random(): [0.0, 1.0)
			int t = (int) Math.floor(Math.random()*Size_Iu_i);
			int i = Iu_i.get(t);                                                       // === i
			// ----------------------------- Iu_j --------------------------------
			List<Integer> Iu_j = new ArrayList<Integer>(Iu_j_data.get(u));			   // === Iu_j
			int Size_Iu_j = Iu_j.size();
			// ----------------------------- Iu_p --------------------------------
			List<Integer> Iu_p = new ArrayList<Integer>(Iu_p_data.get(u));     		   // === Iu_p
			int Size_Iu_p = Iu_p.size();
			// ----------------------------- Iu_c --------------------------------
			List<Integer> Iu_c = new ArrayList<Integer>(Iu_c_data.get(u));           // === Iu_c
			int Size_Iu_c = Iu_c.size();
			// =========================================================
			// === for each user u: the size of p or j may equal to zero, 
			// === so there are 4 situations in which items can be divided. 
			// [1] i,p,c,j; 4 pairwise preferences: Uic, Uij, Upc, Upj
			// [2] i,c,j; 2 pairwise preferences: Uic, Uij
			// [3] i,p,c; 2 pairwise preferences: Uic, Upc
			// [4] i,c; 1 pairwise preferences: Uic
			// =========================================================
			if (Size_Iu_p == 0 && Size_Iu_j == 0){ 
				total_loss += train_ic(Size_Iu_c, Iu_c, u, i);                                        // [4]
			}
			else if (Size_Iu_p != 0 && Size_Iu_j == 0){
				total_loss += train_ipc(Size_Iu_p, Iu_p, Size_Iu_c, Iu_c, u, i);                      // [3]
			}
			else if (Size_Iu_p == 0 && Size_Iu_j != 0){
				total_loss += train_icj(Size_Iu_j, Iu_j, Size_Iu_c, Iu_c, u, i);                      // [2]
			}
			else {
				total_loss += train_ipcj(Size_Iu_j, Iu_j, Size_Iu_p, Iu_p, Size_Iu_c, Iu_c, u, i);    // [1]
			}
		}
		if (iter % 100 == 0) {
			total_loss = -total_loss / (float) n / (float) 100;
			System.out.println("iter " + iter + ", total_loss: " + total_loss);
			pre_total_loss = total_loss;
			total_loss = 0f;
		} 
        }
	}
	// =============================================================
	public static float train_ic(Integer Size_Iu_c, List<Integer> Iu_c, Integer u, Integer i){
		float loss_ic = 0;
		// --- randomly sample an item $c$
		int t4 = (int) Math.floor(Math.random()*Size_Iu_c);
		int c = Iu_c.get(t4);                                                // === c
		// ----------------------------------------------------------------------------
		// --- calculate loss_
		float r_ui = biasV[i];
		float r_uc = biasV[c];
		for (int f=0; f<d; f++){
			r_ui += U[u][f] * V[i][f];			
			r_uc += U[u][f] * V[c][f];	
		}
		float r_Xuic = r_ui - r_uc;
		float loss_uic =  1f / (1f + (float) Math.pow(Math.E, r_Xuic));
		// --------------------------------------------------------------------------
		// --- update $U_{w\cdot}$
		for(int f=0; f<d; f++)
		{
			U[u][f] += gamma * ( loss_uic * (V[i][f] - V[c][f]) - alpha_u * U[u][f] );
		}
		// ---------------------------------------------------
		// --- update $V_{i\cdot}$
		for (int f=0; f<d; f++)
		{
			V[i][f] += gamma * ( loss_uic * U[u][f] - alpha_v * V[i][f] );
		}
		// --- update $V_{c\cdot}$
		for (int f=0; f<d; f++)
		{
			V[c][f] += gamma * ( loss_uic * U[u][f] * (-1) - alpha_v * V[c][f] );
		}
		// --- update $b_i$
		biasV[i] += gamma * ( loss_uic - beta_v * biasV[i] );
		// --- update $b_c$
		biasV[c] += gamma * ( loss_uic * (-1) - beta_v * biasV[c] );
		// ---------------------------------------------------
		float bpr_loss = (float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuic)));
		float reg_loss = (float) (0.5*alpha_u*(norm(U[u])+norm(V[i])+norm(V[c])+biasV[i]*biasV[i]+biasV[c]*biasV[c]));
		loss_ic = bpr_loss - reg_loss;

		return loss_ic;
	}

	public static float train_ipc(Integer Size_Iu_p, List<Integer> Iu_p, Integer Size_Iu_c, List<Integer> Iu_c, Integer u, Integer i){
		float loss_ipc = 0;
		// --- randomly sample an item $l$
		int t3 = (int) Math.floor(Math.random()*Size_Iu_p);
		int p = Iu_p.get(t3);                                                  // === p
		// --- randomly sample an item $c$
		int t4 = (int) Math.floor(Math.random()*Size_Iu_c);
		int c = Iu_c.get(t4);                                                // === c
		// ------------------------------------------------------------------------------------------------
		// --- calculate loss_
		float r_ui = biasV[i];
		float r_up = biasV[p];
		float r_uc = biasV[c];		
		for (int f=0; f<d; f++)
		{
			r_ui += U[u][f] * V[i][f];
			r_up += U[u][f] * V[p][f];
			r_uc += U[u][f] * V[c][f];
		}
		float r_Xuic = r_ui - r_uc;
		float r_Xupc = r_up - r_uc;

		float loss_uic =  1f / (1f + (float) Math.pow(Math.E, r_Xuic) );
		float loss_ulc =  1f / (1f + (float) Math.pow(Math.E, r_Xupc) );
		// --------------------------------------------------------------------------
		// --- update $U_{w\cdot}$
		for(int f=0; f<d; f++)
		{
			U[u][f] += gamma * ( loss_uic * (V[i][f] - V[c][f]) + loss_ulc * (V[p][f] - V[c][f]) - alpha_u * U[u][f] );
		}
		// ---------------------------------------------------
		// --- update $V_{i\cdot}$
		for (int f=0; f<d; f++)
		{
			V[i][f] += gamma * ( loss_uic * U[u][f] - alpha_v * V[i][f] );
		}
		// --- update $V_{p\cdot}$
		for (int f=0; f<d; f++)
		{
			V[p][f] += gamma * ( loss_ulc * U[u][f] - alpha_v * V[p][f] );
		}
		// --- update $V_{c\cdot}$
		for (int f=0; f<d; f++)
		{
			V[c][f] = V[c][f] + gamma * ( (loss_uic + loss_ulc) * U[u][f] * (-1) - alpha_v * V[c][f] );
		}
		// ---------------------------------------------------
		// --- update $b_i$
		biasV[i] += gamma * ( loss_uic - beta_v * biasV[i] );
		// --- update $b_p$
		biasV[p] += gamma * ( loss_ulc - beta_v * biasV[p] );
		// --- update $b_c$
		biasV[c] += gamma * ( (loss_uic + loss_ulc) * (-1) - beta_v * biasV[c] );
		// ---------------------------------------------------
		float bpr_loss = (float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuic)))+(float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xupc)));
		float reg_loss = (float) (0.5*alpha_u*(norm(U[u])+norm(V[i])+norm(V[p])+norm(V[c])+biasV[i]*biasV[i]+biasV[p]*biasV[p]+biasV[c]*biasV[c]));
		loss_ipc = bpr_loss - reg_loss;

		return loss_ipc;
	}

	public static float train_icj(Integer Size_Iu_j, List<Integer> Iu_j, Integer Size_Iu_c, List<Integer> Iu_c, Integer u, Integer i){
		float loss_icj = 0;
		// --- randomly sample an item $j$
		int t1 = (int) Math.floor(Math.random()*Size_Iu_j);
		int j = Iu_j.get(t1);                                                  // === j
		// --- randomly sample an item $c$
		int t2 = (int) Math.floor(Math.random()*Size_Iu_c);
		int c = Iu_c.get(t2);                                                 // === c
		// -------------------------------------------------------------------------------------
		// ---calculate loss_
		float r_ui = biasV[i];
		float r_uc = biasV[c];
		float r_uj = biasV[j];
		for (int f=0; f<d; f++)
		{
			r_ui += U[u][f] * V[i][f];
			r_uc += U[u][f] * V[c][f];
			r_uj += U[u][f] * V[j][f];
		}
		float r_Xuic = r_ui - r_uc;
		float r_Xuij = r_ui - r_uj;				

		float loss_uic =  1f / (1f + (float) Math.pow(Math.E, r_Xuic) );
		float loss_uij =  1f / (1f + (float) Math.pow(Math.E, r_Xuij) );	

		// ---------------------------------------------------------------------
		// --- update $U_{w\cdot}$
		for(int f=0; f<d; f++)
		{
			U[u][f] = U[u][f] + gamma * ( loss_uic * (V[i][f] - V[c][f]) + loss_uij * (V[i][f] - V[j][f]) - alpha_u * U[u][f] );
		}
		// ---------------------------------------------------
		// --- update $V_{i\cdot}$
		for (int f=0; f<d; f++)
		{
			V[i][f] = V[i][f] + gamma * ( (loss_uic + loss_uij) * U[u][f] - alpha_v * V[i][f] );
		}
		// --- update $V_{c\cdot}$
		for (int f=0; f<d; f++)
		{
			V[c][f] = V[c][f] + gamma * ( loss_uic * U[u][f] * (-1) - alpha_v * V[c][f] );
		}
		// --- update $V_{j\cdot}$
		for (int f=0; f<d; f++)
		{
			V[j][f] = V[j][f] + gamma * ( loss_uij * U[u][f] * (-1) - alpha_v * V[j][f] );
		}
		// ---------------------------------------------------
		// --- update $b_i$
		biasV[i] = biasV[i] + gamma * ( loss_uic + loss_uij - beta_v * biasV[i] );
		// --- update $b_c$
		biasV[c] = biasV[c] + gamma * ( loss_uic * (-1) - beta_v * biasV[c] );
		// --- update $b_j$
		biasV[j] = biasV[j] + gamma * ( loss_uij * (-1) - beta_v * biasV[j] );
		// ---------------------------------------------------
		float bpr_loss = (float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuic))) + (float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuij)));
		float reg_loss = (float) (0.5*alpha_u*(norm(U[u])+norm(V[i])+norm(V[c])+norm(V[j])+biasV[i]*biasV[i]+biasV[c]*biasV[c]+biasV[j]*biasV[j]));
		loss_icj = bpr_loss - reg_loss;

		return loss_icj;
	}
	
	public static float train_ipcj(Integer Size_Iu_j, List<Integer> Iu_j, Integer Size_Iu_p, List<Integer> Iu_p, Integer Size_Iu_c, List<Integer> Iu_c, Integer u, Integer i){
		float loss_ipcj = 0;
		// --- randomly sample an item $j$
		int t1 = (int) Math.floor(Math.random()*Size_Iu_j);
		int j = Iu_j.get(t1);                                                  // === j
		// --- randomly sample an item $p$
		int t3 = (int) Math.floor(Math.random()*Size_Iu_p);
		int p = Iu_p.get(t3);                                                  // === p
		// --- randomly sample an item $c$
		int t4 = (int) Math.floor(Math.random()*Size_Iu_c);
		int c = Iu_c.get(t4);                                                // === c
		// ------------------------------------------------------------------------------------------------
		// --- calculate loss_
		float r_ui = biasV[i];
		float r_up = biasV[p];
		float r_uc = biasV[c];
		float r_uj = biasV[j];
		for (int f=0; f<d; f++)
		{
			r_ui += U[u][f] * V[i][f];
			r_up += U[u][f] * V[p][f];
			r_uc += U[u][f] * V[c][f];
			r_uj += U[u][f] * V[j][f];
		}
		float r_Xuic = r_ui - r_uc;
		float r_Xuij = r_ui - r_uj;
		float r_Xupc = r_up - r_uc;
		float r_Xupj = r_up - r_uj;

		float loss_uik =  1f / (1f + (float) Math.pow(Math.E, r_Xuic) );
		float loss_uij =  1f / (1f + (float) Math.pow(Math.E, r_Xuij) );
		float loss_upk =  1f / (1f + (float) Math.pow(Math.E, r_Xupc) );
		float loss_upj =  1f / (1f + (float) Math.pow(Math.E, r_Xupj) );
		// --------------------------------------------------------------------------
		// --- update $U_{w\cdot}$
		for(int f=0; f<d; f++)
		{
			U[u][f] = U[u][f] + gamma * ( loss_uik * (V[i][f] - V[c][f]) + loss_uij * (V[i][f] - V[j][f])  + loss_upk * (V[p][f] - V[c][f]) + loss_upj * (V[p][f] - V[j][f])- alpha_u * U[u][f] );
		}
		// ---------------------------------------------------
		// --- update $V_{i\cdot}$
		for (int f=0; f<d; f++)
		{
			V[i][f] = V[i][f] + gamma * ( (loss_uik + loss_uij) * U[u][f] - alpha_v * V[i][f] );
		}
		// --- update $V_{p\cdot}$
		for (int f=0; f<d; f++)
		{
			V[p][f] = V[p][f] + gamma * ( (loss_upk + loss_upj) * U[u][f] - alpha_v * V[p][f] );
		}
		// --- update $V_{c\cdot}$
		for (int f=0; f<d; f++)
		{
			V[c][f] = V[c][f] + gamma * ( (loss_uik + loss_upk) * U[u][f] * (-1) - alpha_v * V[c][f] );
		}
		// --- update $V_{j\cdot}$
		for (int f=0; f<d; f++)
		{
			V[j][f] = V[j][f] + gamma * ( (loss_uij + loss_upj) * U[u][f] * (-1) - alpha_v * V[j][f] );
		}
		// ---------------------------------------------------
		// --- update $b_i$
		biasV[i] = biasV[i] + gamma * ( loss_uik + loss_uij - beta_v * biasV[i] );
		// --- update $b_p$
		biasV[p] = biasV[p] + gamma * ( loss_upk + loss_upj - beta_v * biasV[p] );
		// --- update $b_c$
		biasV[c] = biasV[c] + gamma * ( (loss_uik + loss_upk) * (-1) - beta_v * biasV[c] );
		// --- update $b_j$
		biasV[j] = biasV[j] + gamma * ( (loss_uij + loss_upj) * (-1) - beta_v * biasV[j] );
		// ---------------------------------------------------
		float bpr_loss = (float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuic)))+(float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xuij)))+(float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xupc)))+(float) Math.log(1f / (1f + (float) Math.pow(Math.E, -r_Xupj)));
		float reg_loss = (float) (0.5*alpha_u*(norm(U[u])+norm(V[i])+norm(V[p])+norm(V[c])+norm(V[j])+biasV[i]*biasV[i]+biasV[p]*biasV[p]+biasV[c]*biasV[c]+biasV[j]*biasV[j]));
		loss_ipcj = bpr_loss - reg_loss;
		// ===================================================
		return loss_ipcj;
	}
	
	// =============================================================
	public static double consine(Collection<Integer> list1, Collection<Integer> list2) {
		
		double result = 0;
		int num;
		int sum = 0;
		List<Integer> line1 = new ArrayList<Integer>(list1);
		List<Integer> line2 = new ArrayList<Integer>(list2);
		
		num = list1.size();
		for(int i=0; i<num ; i++)
		{
			sum += line1.get(i)*line2.get(i); 
		}
		
		double a,b;
		a = sum(line1);
		b = sum(line2);
		a = Math.sqrt(a);
		b = Math.sqrt(b);
		result = a*b;
		if(result == 0) 
		{
			result = 0;
		}
		else 
		{
			result = sum/result;
		}
		return result;
	}
    
    public static double sum(Collection<Integer> list1) {
		
		double result = 0;
		int num;
		
		
		List<Integer> line1 = new ArrayList<Integer>(list1);
		num = list1.size();
		for(int i=0; i<num; i++)
		{
			result += line1.get(i)*line1.get(i);
		}

		return result;	
	}
	public static float norm(float[] aa) {
    	float res = 0f;
    	for(float i :aa) {
    		res += i * i;
    	}
    	return (float) Math.sqrt(res);
    }
    
    // =============================================================
    public static void testRanking(HashMap<Integer, HashSet<Integer>> TestData)
    {
		// TestData: user->items 
		// ==========================================================
		float[] PrecisionSum = new float[topK+1];
		float[] RecallSum = new float[topK+1];	
		float[] F1Sum = new float[topK+1];
		float[] NDCGSum = new float[topK+1];
		float[] OneCallSum = new float[topK+1];
		float MRRSum = 0;
		float MAPSum = 0;
		// --- calculate the best DCG, which can be used later
		float[] DCGbest = new float[topK+1];
		for (int k=1; k<=topK; k++)
		{
			DCGbest[k] = DCGbest[k-1];
			DCGbest[k] += 1/Math.log(k+1); 
		}
		
		// --- number of test cases
    	int UserNum_TestData = TestData.keySet().size(); 
    	
    	for(int u=1; u<=n; u++)
    	{
    		// --- check whether the user $u$ is in the test set
    		if (!TestData.containsKey(u))
    			continue;
    		
    		// ---
    		Set<Integer> ItemSet_u_TrainData = new HashSet<Integer>(); 
    		if (TrainData.containsKey(u))
    		{
    			ItemSet_u_TrainData = TrainData.get(u).keySet();
    		}
    		HashSet<Integer> ItemSet_u_TestData = TestData.get(u); 

    		// --- the number of preferred items of user $u$ in the test data 
    		int ItemNum_u_TestData = ItemSet_u_TestData.size();   
    		
    		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    		
    		// --- prediction    		
    		HashMap<Integer, Float> item2Prediction = new HashMap<Integer, Float>();
    		item2Prediction.clear();
    		
    		for(int i=1; i<=m; i++)
    		{
    			// --- (1) check whether item $i$ is in the whole item set
    			// --- (2) check whether item $i$ appears in the training set of user $u$
    			// --- (3) check whether item $i$ is in the ignored set of items    			
    			if ( !ItemSetWhole.contains(i) || ItemSet_u_TrainData.contains(i) )  
    				continue;

    			// --- prediction via inner product
        		float pred = biasV[i];
        		for (int f=0; f<d; f++)
        		{
        			pred += U[u][f]*V[i][f];
        		}       			
        		item2Prediction.put(i, pred);
        	}
    		// --- sort
    		List<Map.Entry<Integer,Float>> listY = 
    				new ArrayList<Map.Entry<Integer,Float>>(item2Prediction.entrySet()); 		
    		Collections.sort(listY, new Comparator<Map.Entry<Integer,Float>>()
    		{
    			public int compare( Map.Entry<Integer, Float> o1, Map.Entry<Integer, Float> o2 )   
    			{
    				return o2.getValue().compareTo( o1.getValue() ); 
    			}
    		});
    		
    		// ===========================================================
    		// === Evaluation: TopK Result 
    		// --- Extract the topK recommended items
    		int k=1;
    		int[] TopKResult = new int [topK+1];    		
    		Iterator<Entry<Integer, Float>> iter = listY.iterator();
    		while (iter.hasNext())
    		{
    			if(k>topK)
    				break;
    			
    			Map.Entry<Integer, Float> entry = (Map.Entry<Integer, Float>) iter.next(); 
    			int itemID = entry.getKey();
    			TopKResult[k] = itemID;
    			k++;
    		}
    		// --- TopK evaluation
    		int HitSum = 0;
    		float[] DCG = new float[topK+1];
    		float[] DCGbest2 = new float[topK+1];
    		for(k=1; k<=topK; k++)
    		{
    			// ---
    			DCG[k] = DCG[k-1];
    			int itemID = TopKResult[k];
    			if ( ItemSet_u_TestData.contains(itemID) )
    			{
        			HitSum += 1;
        			DCG[k] += 1 / Math.log(k+1);
    			}
    			// --- precision, recall, F1, 1-call
    			float prec = (float) HitSum / k;
    			float rec = (float) HitSum / ItemNum_u_TestData;    			
    			float F1 = 0;
    			if (prec+rec>0)
    				F1 = 2 * prec*rec / (prec+rec);
    			PrecisionSum[k] += prec;
    			RecallSum[k] += rec;
    			F1Sum[k] += F1;
    			// --- in case the the number relevant items is smaller than k 
    			if (ItemSet_u_TestData.size()>=k)
    				DCGbest2[k] = DCGbest[k];
    			else
    				DCGbest2[k] = DCGbest2[k-1];
    			NDCGSum[k] += DCG[k]/DCGbest2[k];
    			// ---
    			OneCallSum[k] += HitSum>0 ? 1:0; 
    		}
    		// ===========================================================
    		
    		// ===========================================================
    		// === Evaluation: Reciprocal Rank
    		if (flagMRR)
    		{
	    		int p = 1;
	    		iter = listY.iterator();    		
	    		while (iter.hasNext())
	    		{	
	    			Map.Entry<Integer, Float> entry = (Map.Entry<Integer, Float>) iter.next(); 
	    			int itemID = entry.getKey();
	    			
	    			// --- we only need the position of the first relevant item
	    			if(ItemSet_u_TestData.contains(itemID))    				
	    				break;
	
	    			p += 1;
	    		}
	    		MRRSum += 1 / (float) p;
    		}
    		// ===========================================================
    		
    		// ===========================================================
    		// === Evaluation: Average Precision
    		if (flagMAP)
    		{
	    		int p = 1; // the current position
	    		float AP = 0;
	    		int HitBefore = 0; // number of relevant items before the current item
	    		iter = listY.iterator();    		
	    		while (iter.hasNext())
	    		{	
	    			Map.Entry<Integer, Float> entry = (Map.Entry<Integer, Float>) iter.next(); 
	    			int itemID = entry.getKey();
	    			
	    			if(ItemSet_u_TestData.contains(itemID))
	    			{
	    				AP += 1 / (float) p * (HitBefore + 1);
	    				HitBefore += 1;
	    			}
	    			p += 1;
	    		}
	    		MAPSum += AP / ItemNum_u_TestData;
    		}
    	}
    	
    	// =========================================================
    	// --- precision@k
    	for(int k=topK; k<=topK; k++)
    	{
    		float prec = PrecisionSum[k]/UserNum_TestData;
    		Pre_ave += prec;
    		System.out.println("Prec@"+Integer.toString(k)+":"+Float.toString(prec));   
    	}
    	// --- recall@k
    	for(int k=topK; k<=topK; k++)
    	{
    		float rec = RecallSum[k]/UserNum_TestData;
    		Rec_ave += rec;
    		System.out.println("Rec@"+Integer.toString(k)+":"+Float.toString(rec));  
    	}
    	// --- F1@k
    	for(int k=topK; k<=topK; k++)
    	{
    		float F1 = F1Sum[k]/UserNum_TestData;
    		F1_ave += F1;
    		System.out.println("F1@"+Integer.toString(k)+":"+Float.toString(F1));   
    	}
    	// --- NDCG@k
    	for(int k=topK; k<=topK; k++)
    	{
    		float NDCG = NDCGSum[k]/UserNum_TestData;
    		NDCG_ave += NDCG;
    		System.out.println("NDCG@"+Integer.toString(k)+":"+Float.toString(NDCG));
    	}
    	// --- MRR
    	float MRR = MRRSum/UserNum_TestData;
    	MRR_ave += MRR;
    	System.out.println("MRR:" + Float.toString(MRR));
    	// --- MAP
    	float MAP = MAPSum/UserNum_TestData;
    	MAP_ave += MAP;
    	System.out.println("MAP:" + Float.toString(MAP));
    }
    
}
	
