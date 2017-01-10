import java.util.*;
import java.io.*;

public class CorrelationClient {
 private static double pVal, tau;

    private static void handleNewPriority() {
        Scanner in = new Scanner(System.in);
        System.out.println("Computations for gene of interest has been completed... Enter new gene of interest if desired: ");
        try {
            String newPriority = in.next();
            if (!Correlate.hasCompleted()) {
                Correlate.setPriorityGene(newPriority);
            }
        }catch (Exception e) {
            ///
        }
    }
    
    public static void main(String[] args) {
        final Scanner in = new Scanner(System.in);
        Correlate.Listener listener = new Correlate.Listener() {
            public void onPriorityCompleted(Correlate.EdgeList edgeList) {
                System.out.println("--------------------------------------------------");
                for (int i = 0; i < edgeList.first().length; i++) {
                    String first = edgeList.getFirstName(i);
                    String second = edgeList.getSecondName(i);
                    float value = edgeList.values()[i];
                    System.out.format("%-10s%-10s%-10.2f\n", first, second, value);
                }
                if (!Correlate.hasCompleted()) {
                  Correlate.graphData(edgeList);
                    handleNewPriority();
                }
            }
            
            public void onCompleted() {
                try {
                    Correlate.EdgeList list = Correlate.getEdgeList();
                    Correlate.graphData(list);
                } catch (Exception e) {
                    //
                }
                System.out.println("All computations have completed");
                while (true) {
                    System.out.println("Enter absolute file path to save correlations: ");
                    String corrFile = in.nextLine();
                    try {
                        boolean corrSave = Correlate.outCorrMatrix(corrFile);
                        if (corrSave) break;
                        else System.out.println("Error in saving correlation...");
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                }
                while (true) {
                    System.out.println("Enter absolute file path to save p-values: "); 
                    String pvalFile = in.nextLine();
                    try {
                        boolean pvalSave = Correlate.outPvalMatrix(pvalFile);
                        if (pvalSave) break;
                        else System.out.println("Error in saving p-values");
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                }
                while (true) {
                    System.out.println("Enter absolute file path to save edge list: ");
                    String edgeFile = in.nextLine();
                    try {
                        boolean edgeSave = Correlate.outEdgeList(edgeFile);
                        if (edgeSave) break;
                        else System.out.println("Error in saving edge list");
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                }
                System.out.println("Results have been successfully saved");
                Correlate.end();
            }
            
            public void onFailed() {
                System.out.println("Computations have failed");
                Correlate.end();
            }
        };
        
        System.out.println("Please input the filepath of your aggregated data set:");
        String filePath = in.nextLine();
        
        if (Correlate.isInitialized()) {
            Correlate.end();
        }
        
        File file = new File(filePath);
        try {
            Correlate.start(file);
        } catch (Exception e) {
            System.out.println(e);
            return;
        }
        
        System.out.print("Enter your tau threshold: ");
        tau = Double.parseDouble(in.nextLine().trim());
        
        System.out.print("Enter your significance level: ");
        pVal = Double.parseDouble(in.nextLine().trim());
        
        String priorityGene = null;
        while (true) {
            System.out.print("Do you have a specific gene of interest? (Y/N): ");
            String val = in.nextLine().trim();
            if (val.equals("Y")) {
                System.out.print("Enter your gene of interest: ");
                priorityGene = in.nextLine();
                break;
            } else if (val.equals("N")) {
                break;
            } else {
                System.out.println("Error: Incorrect input - please enter Y or N");
            }
        }
        
        try {
            Correlate.corrData(pVal, tau, priorityGene, listener);
        } catch (Exception e) {
            System.out.println(e);
            return;
        }
        
        System.out.println("The correlational computations are running...");
    }
}