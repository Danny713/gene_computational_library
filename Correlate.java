/*import libraries*/
import java.io.*;
import java.util.*;
import org.rosuda.JRI.Rengine;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RBool;
import javax.swing.*;

/* Some of the logic has been taken from the GenEx Project */

/* Static class to conduct correlational computations */
public class Correlate {
    private static final String SAFETY = "A";
    
    private static Rengine re;
    private static boolean isInitialized;
    
    private static boolean isAggregated;
    private static boolean isCorrelated;
    
    private static String aggDT;
    private static String corrEnv;
    
    private static String adjEnv;
    
    private static EdgeList edgeList;
    private static String[] origGeneList;
    
    private static File filePath;
    private static Listener listener;
    private static String priorityGene;
    private static double pVal;
    private static double tau;
    
    private static Thread prevThread;
    
    /* -------------------------- Methods to interface with R engine -------------------------- */
    /*start Rengine*/
    private static Rengine startR() {
        Rengine re = Rengine.getMainEngine();
        //check if R Engine exists
        if (re == null) {
            // Create R Engine
            String[] Rargs = {"--vanilla"};
            re = new Rengine(Rargs, false, null);
            
            re.eval("library(data.table)");
            re.eval("library(WGCNA)");
            re.eval("library(rJava)");

            re.eval("source('./Backend/R/queue.R')");
            re.eval("source('./Backend/R/corrData.R')");
            re.eval("dyn.load('./Backend/src/pvalAndTriFunctions.so')");
            
            re.eval("gc()");
        }
        
        return re;
    }
    
    /*ends R session -- only include if you're done with entire Java session*/
    private static void endR(Rengine re) {
        if (re != null) {
            re.end();
        }
    }
    
    /*check if R object is NOT null*/
    private static boolean notNull(String object, Rengine re) {
        if (re == null)
            return false;
        //if not null, re.eval(check).asBool() == FALSE; .isFALSE() == TRUE
        return re.eval(String.format("is.null(%s)", object)).asBool().isFALSE();
    }
    
    /* ----------------------------------------------------------------------------------------- */
    private static String getExtension(File file) {
        return getExtension(file.getName());
    }
    
    private static String getExtension(String filename) {
        String ext = null;
        int i = filename.lastIndexOf('.');
        if (i > 0 &&  i < filename.length() - 1) {
            ext = filename.substring(i).toLowerCase();
        }
        return ext;
    }
    
    private static String createFileName(String filename) {
        int i = filename.lastIndexOf('.');
        String name;
        if (i > 0 &&  i < filename.length() - 1) {
            name = filename.substring(0, i).toLowerCase();
        } else {
            name = "temporary";
        }
        return SAFETY + name;
    }
    
    private static boolean loadAggDT() {
        String ext = getExtension(filePath);
        
        // only handle .csv files
        if (!ext.equals(".csv")) {
            return false;
        }   

        String aggregate = String.format("%s <- loadAggDT(\'%s\')", aggDT, filePath.getAbsolutePath());
        re.eval(aggregate);

        isAggregated = notNull(aggDT, re);
        return isAggregated;
    }

    /* Initiates the interface using the aggregated data set, filePath. */
    public static boolean start(File filePath) throws Exception {
        if (isInitialized) {
            throw new Exception("Correlational computations are already running");
        }
        re = startR();
        if (re == null) {
            end();
            return false;
        }

        isInitialized = true;
        Correlate.filePath = filePath;
        
        aggDT = createFileName(filePath.getName()) + "_aggDT";
        if (!loadAggDT()) {
            end();
            return false;
        }
        return true;
    }

    /* Determines if a session is running */
    public static boolean isInitialized() {
        return isInitialized;
    }
    
    /* Determines if all computations have completed */
    public static boolean hasCompleted() {
        return isCorrelated;
    }
    
    /* Ends the entire session */
    public static void end() {
        if (isInitialized) {
            endR(re);
            isInitialized = false;
            isAggregated = false;
            isCorrelated = false;
            aggDT = null;
            
            adjEnv = null;
            
            corrEnv = null;
            edgeList = null;
            origGeneList = null;
            filePath = null;
            listener = null;
            priorityGene = null;
            prevThread = null;
            pVal = -1;
            tau = -1;
        }
    }

    /* Initiates computations based on correlation coefficient threshold, tau,
     * significance value, pval. Listener is used to notify client */
    public static void corrData(final double pval, final double tau, final Listener listener) throws Exception {
        corrData(pval, tau, null, listener);
    }
    
    /* Initiate computation with an initial gene of interest, priorityGene */
    public static void corrData(final double pval, final double tau, String priorityGene, final Listener listener) throws Exception {
        if (isAggregated) {
            Correlate.priorityGene = priorityGene;
            Correlate.pVal = pval;
            Correlate.tau = tau;
            
            assert listener != null;
            Correlate.listener = listener;
            // call corrData in a background thread
            Runnable task = new Runnable() {
                @Override
                public void run() {
                    String correlate = String.format("%s <- corrData(%s, %f, %f)", corrEnv, aggDT, pval, tau);
                    re.eval(correlate);
                    //check if successful
                    isCorrelated = notNull(corrEnv, re);
                    if (prevThread != null && prevThread.isAlive()) {
                        prevThread.interrupt();
                    }
                    
                    if (isCorrelated) {
                        listener.onCompleted();
                    } else {
                        listener.onFailed();
                    }
                }
            };
            
            Thread thread = new Thread(task);
            thread.start();
        } else {
            end();
            throw new Exception("No correlational computations are running");
        }
    }
    
    /* Set the gene of interest to priority */
    public static void setPriorityGene(String priority) throws Exception {
        if (!isInitialized) throw new Exception("No correlational computations are running");
        priorityGene = priority;
    }
    
    /* Returns all co-expressed genes with their corresponding p-values. This method should only
     * be called once all computations have concluded. */
    public static EdgeList getEdgeList() throws Exception {
        if (!isInitialized) 
            throw new Exception("No correlational computations are running");
        if (!isCorrelated)
            throw new Exception("Correlational computations have failed or are still running");
        
        return edgeList;
    }
    
    /* Returns all pairs of co-expressed genes, where at least one of the two genes are in the
     * geneList. This method can be called before all computations have concluded. */
    public static EdgeList getEdgeList(String[] geneList) throws Exception {
        if (!isInitialized) {
            throw new Exception("No correlational computations are running");
        }
        HashSet<String> geneSet = new HashSet<String>(Arrays.asList(geneList));
        ArrayList<Integer> indexList = new ArrayList<Integer>();
        for (int i = 0; i < edgeList.first().length; i++) {
            if (geneSet.contains(edgeList.getFirstName(i)) || geneSet.contains(edgeList.getSecondName(i))) {
                indexList.add(i);
            }
        }
        int[] first = new int[indexList.size()];
        int[] second = new int[indexList.size()];
        float[] values = new float[indexList.size()];
        
        int arrayIndex = 0;
        for (Integer index : indexList) {
            first[arrayIndex] = edgeList.first[index];
            second[arrayIndex] = edgeList.second[index];
            values[arrayIndex] = edgeList.values[index];
            arrayIndex++;
        }
        
        return new EdgeList(first, second, edgeList.names, values);
    }
    
    /* Construct a gene co-expression network using the provided edgelist */
    public static void graphData(EdgeList edgeList) {
        try {
            SingleGraph sg = new SingleGraph(edgeList, "hello");
            GraphFramer f = new GraphFramer(sg, pVal, tau);
        } catch (NullPointerException e) {
            System.out.println(e);
        }
    }
    
    /* Save the correlation coefficients to file at filePath */
    public static boolean outCorrMatrix(String filePath) throws Exception {
        boolean success = false;
        
        if (!getExtension(filePath).equals(".txt")) {
            throw new Exception("File is not a txt file");
        }
        if (!isInitialized) {
            throw new Exception("No correlational computations are running");
        }
        
        if (isCorrelated) {
            //output correlation matrix
            String write = String.format("OC <- outMatrix(%s$corrVec, %s$index, \'%s\')", corrEnv, corrEnv, filePath);
            re.eval(write);
            //check if successful
            success = notNull("OC", re);
        }
        return success;
        
    }
    
    /* Save the p-values to file at filePath */
    public static boolean outPvalMatrix(String filePath) throws Exception {
        boolean success = false;
        
        if (!getExtension(filePath).equals(".txt")) {
            throw new Exception("File is not a txt file");
        }
        
        if (!isInitialized) {
            throw new Exception("No correlational computations are running");
        }
        
        if (isCorrelated) {
            String write = String.format("OP <- outMatrix(%s$pVec, %s$index, \'%s\')", corrEnv, corrEnv, filePath);
            re.eval(write);
            //check if successful
            success = notNull("OP", re);
        }
        return success;
    }
    
    /* Save the edge list to file at filePath */
    public static boolean outEdgeList(String filePath) throws Exception {
        boolean success = false;
        
        try {
            File file = new File(filePath);
            PrintWriter writer = new PrintWriter(file);
            EdgeList edgeList = getEdgeList();
            for (int i = 0; i < edgeList.first().length; i++) {
                String first = edgeList.getFirstName(i);
                String second = edgeList.getSecondName(i);
                float value = edgeList.values()[i];
                writer.println(first + "," + second + "," + value);
            }
            writer.close();
            success = true;
        } catch (Exception e) {
            System.out.println(e);
        }

        return success;
    } 
    
    
    /* -------------------------- Methods that only should be called by R -------------------------- */
    
    /* Called by R to retrieve most up-to-date priority gene */
    public static String getPriorityGene() {
        return priorityGene;
    }
    
    /* Should only be accessed by R-backend */
    public static void handleUpdatedEdgeList(Object df[], String names[], boolean isLast) {
        double[] first = (double[]) df[0];
        double[] second = (double[]) df[1];
        double[] values = (double[]) df[2];
        
        int[] newFirst = new int[first.length];
        int[] newSecond = new int[second.length];
        float[] newValues = new float[values.length];
        
        for (int i = 0; i < first.length; i++) {
            newFirst[i] = (int) first[i];
            newSecond[i] = (int) second[i];
            newValues[i] = (float) values[i];
        }
        
        if (origGeneList == null) {
            origGeneList = new String[names.length];
            for (int i = 0; i < names.length; i++) {
                origGeneList[i] = (String) names[i];
            }
        }
        
        edgeList = new EdgeList(newFirst, newSecond, origGeneList, newValues);

        if (!isLast) {
            Runnable task = new Runnable() {
                @Override
                public void run() {
                    listener.onPriorityCompleted(edgeList);
                }
            };
            if (prevThread != null && prevThread.isAlive()) {
                prevThread.interrupt();
            }
            prevThread = new Thread(task);
            prevThread.start();
        }
    }
    
    /* ----------------------------------------------------------------------------------------- */
    
    /**
     * Nested class: EdgeList
     * Contains edge list information
     **/
    public static class EdgeList {
        /*instance variables*/
        private int[] first;
        private int[] second;
        private String[] names;
        private float[] values;
        private int E; //number of edges
        private boolean isValid = false;
        
        /*constructor*/
        public EdgeList(int[] first, int[] second, String[] names, float[] values) {
            this.first = first;
            this.second = second;
            this.names = names;
            this.values = values;
            
            E = values.length;
            
            //make sure parallel arrays (first, second, values) are same size
            if ((E == first.length) && (E == second.length)) {
                isValid = true;
            }
        }
        
        /*check if valid*/
        public boolean isValid() {
            return isValid;
        }
        
        /*access methods*/
        public int[] first() {
            return first;
        }
        public String getFirstName(int i) {
            return names[first[i] - 1];
            //this may throw array out of bounds error
            //subtract 1 because R is 1-indexed and Java is 0-indexed
        }
        
        public int[] second() {
            return second;
        } 
        public String getSecondName(int i) {
            return names[second[i] - 1];
            //this may throw array out of bounds error
            //subtract 1 because R is 1-indexed and Java is 0-indexed
        }
        
        public String[] names() {
            return names;
        }
        public float[] values() {
            return values;
        } 
        
        /*number of edges*/
        public int E() {
            return E;
        }

        public List<String> getAllEntries() {
            ArrayList<String> list = new ArrayList<String>();
            for (int i = 0; i < edgeList.first().length; i++) {
                String first = edgeList.getFirstName(i);
                String second = edgeList.getSecondName(i);
                float value = edgeList.values()[i];
                list.add(first + "," + second + "," + value);
            }

            Collections.sort(list);
            return list;
        }
    }
    
    public static interface Listener {
        public void onPriorityCompleted(EdgeList edgelist);
        public void onCompleted();
        public void onFailed();
    }
}

