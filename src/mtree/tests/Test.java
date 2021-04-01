package mtree.tests;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;

import mtree.utils.FibonacciHeap;
import outlierdetection.OD_Query;

public class Test {
    
    public static ArrayList<OD_Query> all_queries = new ArrayList<>();
    
    public static ArrayList<FibonacciHeap.Node<Integer>> l = new ArrayList<>();
    
    public static ArrayList<Integer> getSlides(int currentTime, int w_start) {
        ArrayList<Integer> allOffsets = new ArrayList<>();
        allOffsets.add(w_start);
        for (OD_Query q : all_queries) {
            int n = 1;
            while (n * q.S - q.W < 0) {
                int offset = currentTime + q.S * n - q.W;
                if (!allOffsets.contains(offset)) {
                    allOffsets.add(offset);
                }
                n++;
            }
        }
        Collections.sort(allOffsets);
        return allOffsets;
        
    }

    public static void main(String[] args) {

        //tao
        double[] r_list = new double[]{0.19, 0.95, 1.9, 9.5, 19};
        int[] k_list = new int[]{10, 50, 100, 200, 500};
        int[] w_list = new int[]{700,2000};
        int[] s_list = new int[]{100,  250};
        for (Double r : r_list) {
            for (Integer k : k_list) {
                for (Integer w : w_list) {
                    for (Integer slide : s_list) {
                        all_queries.add(new OD_Query(r, k, w, slide));
                        
                    }
                }
            }
        }
        ArrayList<Integer> slideOffsets = getSlides(4000, 2000);
        for (Integer offset : slideOffsets) {
            System.out.println(offset);
        }
        
    }
    
    public static void readFromFile(String filename) {
        try {
            FileReader fr = new FileReader(new File(filename));
            BufferedReader bfr = new BufferedReader(fr);
            String line = "";
            while ((line = bfr.readLine()) != null) {
                String[] numbers = line.split(" ");
                for (int i = 0; i < 3; i++) {
                    Double number0 = Double.valueOf(numbers[0]);
                    Double number1 = Double.valueOf(numbers[1]);
                    Double number2 = Double.valueOf(numbers[2]);
                }
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Test.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Test.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public static void printHeaps() {
        System.out.println("------------------------------------");
        for (FibonacciHeap.Node<Integer> node : l) {
            
            System.out.print("key=" + node.getKey() + ";next=" + node.next.getKey() + ";prev=" + node.prev.getKey());
            if (node.parent != null) {
                System.out.print(";parent=" + node.parent.getKey());
            } else {
                System.out.print(";parent=" + node.parent);
            }
            if (node.child != null) {
                System.out.println(";child=" + node.child.getKey());
            } else {
                System.out.println(";child=" + node.child);
            }
        }
    }
}
