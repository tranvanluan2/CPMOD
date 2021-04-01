/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mtree.tests;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import mtree.utils.Constants;
import mtree.utils.Utils;
import outlierdetection.CPOD_MQ_Naive;
import outlierdetection.CPOD_MQ_Naive_CompareRK;
import outlierdetection.CPOD_MQ_ShareCore;
import outlierdetection.CPOD_ShareCore_10;
import outlierdetection.CPOD_ShareCore_4;
import outlierdetection.CPOD_ShareCore_5;
import outlierdetection.CPOD_ShareCore_6;
import outlierdetection.CPOD_ShareCore_7;
import outlierdetection.CPOD_ShareCore_8;
import outlierdetection.CPOD_ShareCore_9;
import outlierdetection.MCSKy;
import outlierdetection.New_MQ_CPOD;
import outlierdetection.New_MQ_CPOD_2;
import outlierdetection.OD_Query;
import outlierdetection.SOP;

/**
 *
 * @author luan
 */
public class OD_Test {

    public static int currentTime = 0;
    public static boolean stop = false;
    public static HashSet<Integer> idOutliers = new HashSet<>();
    public static String algorithm;
    public static int numberWindows = 0;
    public static double start;
    public static int w_min = 1000;
    public static int s_min = 50;
    
    public static double overall_avg_outlier_rate=0;
    public static void main(String[] args) {

        readArguments(args);
        for (String arg : args) {
            System.out.print(arg + " ");
        }
        System.out.println("");

        MesureMemoryThread mesureThread = new MesureMemoryThread();
        mesureThread.start();
        Stream s = Stream.getInstance("");

        CPOD_MQ_Naive cmqn = new CPOD_MQ_Naive();
        CPOD_MQ_ShareCore cmsc = new CPOD_MQ_ShareCore();
        CPOD_MQ_Naive_CompareRK cmqnRK = new CPOD_MQ_Naive_CompareRK();
        SOP sop = new SOP();
        New_MQ_CPOD new_cpod = new New_MQ_CPOD();
        CPOD_ShareCore_4 new_cpod4 = new CPOD_ShareCore_4();
        CPOD_ShareCore_5 new_cpod5 = new CPOD_ShareCore_5();
        CPOD_ShareCore_6 new_cpod6 = new CPOD_ShareCore_6();
        CPOD_ShareCore_7 new_cpod7 = new CPOD_ShareCore_7();
        CPOD_ShareCore_8 new_cpod8 = new CPOD_ShareCore_8();
        CPOD_ShareCore_9 new_cpod9 = new CPOD_ShareCore_9();
        CPOD_ShareCore_10 new_cpod10 = new CPOD_ShareCore_10();
        MCSKy mcsky = new MCSKy();
        //create a set of queries
        //fc
//        double[] r_list = new double[]{52.5, 262.5, 525, 2625, 5250};
//        int[] k_list = new int[]{10, 50, 100, 200,  500};
//        int[] k_list = new int[]{10, 30, 50, 70,  100};

//        double[] r_list = new double[]{525};
//        int[] k_list = new int[]{50};
//      
        int seed = 24;
        Random random = new Random();
        random.setSeed(seed);
//        int numR = 50;
//
//        double[] r_list = new double[numR];
//        for (int i = 0; i < numR; i++) {
//            // r_list[i] = 525/2 + 525*1.5*random.nextDouble();//fc
//            r_list[i] = 1.9/2 + 1.9*1.5*random.nextDouble(); //tao
//            //           r_list[i] = 6.5/2 + 6.5*1.5*random.nextDouble(); //em
//            //          r_list[i] = 2.75/2 + 2.75*1.5*random.nextDouble(); //gas
//            //   r_list[i] = 0.45/2 + 0.45*1.5*random.nextDouble(); //stk
////            r_list[i] = 0.028 / 2 + 0.028 * 1.5 * random.nextDouble(); //gau
//        }
        //tao
    //    Double[] r_list = new Double[]{1.9};
        //fc
//        Double[] r_list = new Double[]{525.0};
//        Double[] r_list = new Double[]{6.5};
//        Double[] r_list = new Double[]{0.45};
//        Double[] r_list = new Double[]{0.028};
//        Double[] r_list = new Double[]{2.75};
//        Double[] r_list = new Double[]{115.0};
 //       Integer[] k_list = new Integer[]{50};
//        Integer[] w_list = new Integer[]{10000};
        int num_W = 50;
//        Integer[] w_list = new Integer[num_W];
//        for (int i = 0; i < w_list.length; i++) {
//            w_list[i] = 500 + random.nextInt(60)*500;
//        }
//        int numS = 50;
//        Integer[] s_list = new Integer[numS];
//        for (int i = 0; i < s_list.length; i++) {
//            s_list[i] =  random.nextInt(50)*200;
//        }
//        Integer[] s_list = new Integer[]{500};

//        int[] k_list = new int[]{10, 30, 50, 70,  100};
        //        double[] r_list = new double[]{1.9};
        //        int[] k_list = new int[]{50};
        //        
        //gau
        //        double[] r_list = new double[]{0.0028, 0.014, 0.028, 0.14, 0.28};
        //        int[] k_list = new int[]{10, 50, 100, 200,  500};
        //        int[] k_list = new int[]{10, 30, 50, 70,  100};
        
        ArrayList<Integer> k_list = new ArrayList<>();
        for(int i = 0; i < num_W; i++){
            k_list.add(50/2 + (int)(50*1.5*random.nextDouble()));
        }
        
        ArrayList<Double> r_list = new ArrayList<>();
        for(int i = 0; i < num_W; i++){
            r_list.add(0.028/2 + 0.028*1.5*random.nextDouble());
        }
        
//        Collections.sort(r_list);
//        Collections.sort(k_list);
//        Collections.sort(k_list, Collections.reverseOrder());
        for(int i =0 ; i < num_W; i++){
//            double r = 1.9;
//            double r = 1.9/2 + 1.9*1.5*random.nextDouble();
//            double r = 525.0/2 + 525.0*1.5*random.nextDouble();
//                double r = 525;
//             double r = 115.0;
//            double r = 115.0/2 + 115.0*1.5*random.nextDouble();
//            double r = 6.5/2 + 6.5*1.5*random.nextDouble();
//            double r = 2.75/2 + 2.75*1.5*random.nextDouble();
//            double r = 0.45/2 + 0.45*1.5*random.nextDouble();
//            double r = 0.45/2 + 0.45*1.5*random.nextDouble();
//            double r = 0.028/2 + 0.028*1.5*random.nextDouble();
//            double r = 6.5;
//            double r = 2.75;
//            double r = 0.45;
//            double r = 0.028;
               double r = r_list.get(i);
               int k = k_list.get(i);
//            int k = 50/2 + (int)(50*1.5*random.nextDouble());
//            int k = 50; 
              
//            int k = 50/2 + (int)(50*9*random.nextDouble());
          //  int k = 50/2 + (int)(50*4*random.nextDouble());
          
//            int w = 200 + random.nextInt(150) * 200;
            int w = 100000;
//            int slide = 1600 + random.nextInt(150)*200;
            int slide = 5000;
//            int w = 2000 + random.nextInt(150)*2000;
//            int slide = 16000 + random.nextInt(150)*2000;
            cmqn.add_query(new OD_Query(r, k, w, slide));
            cmsc.add_query(new OD_Query(r, k, w, slide));
            cmqnRK.add_query(new OD_Query(r, k, w, slide));
            sop.add_query(new OD_Query(r, k, w, slide));
            new_cpod4.add_query(new OD_Query(r, k, w, slide));
            new_cpod5.add_query(new OD_Query(r, k, w, slide));
            new_cpod6.add_query(new OD_Query(r, k, w, slide));
            new_cpod7.add_query(new OD_Query(r, k, w, slide));
            new_cpod8.add_query(new OD_Query(r, k, w, slide));
            new_cpod9.add_query(new OD_Query(r, k, w, slide));
            new_cpod10.add_query(new OD_Query(r, k, w, slide));
            mcsky.add_query(new OD_Query(r, k, w, slide));
        }
//        for (Double r : r_list) {
//            for (Integer k : k_list) {
//                for (Integer w : w_list) {
//                    for (Integer slide : s_list) {
//                        cmqn.add_query(new OD_Query(r, k, w, slide));
//                        cmsc.add_query(new OD_Query(r, k, w, slide));
//                        cmqnRK.add_query(new OD_Query(r, k, w, slide));
//                        sop.add_query(new OD_Query(r, k, w, slide));
////                        new_cpod.add_query(new OD_Query(r, k, w, slide));
//                        new_cpod4.add_query(new OD_Query(r, k, w, slide));
//                        new_cpod5.add_query(new OD_Query(r, k, w, slide));
//                        mcsky.add_query(new OD_Query(r, k, w, slide));
//                    }
//                }
//            }
//        }

//        cmqn.add_query(new OD_Query(500, 45));
        double totalTime = 0;
        while (!stop) {
            if (Constants.numberWindow != -1 && numberWindows > Constants.numberWindow) {
                break;
            }
            numberWindows++;
            System.out.println("Num window = " + numberWindows);

            ArrayList<Data> incomingData;
            if (currentTime != 0) {
                incomingData = s.getIncomingData(currentTime, Constants.slide, Constants.dataFile, Constants.matrixType);
                currentTime = currentTime + Constants.slide;
//                System.out.println("Last idx time = " + (incomingData.get(incomingData.size()-1).arrivalTime-1));
            } else {
                incomingData = s.getIncomingData(currentTime, Constants.W, Constants.dataFile, Constants.matrixType);
                currentTime = currentTime + Constants.W;

            }

            start = Utils.getCPUTime(); // requires java 1.5

            /**
             * do algorithm
             */
//            HashMap<OD_Query, HashSet<Data>> result = cmqn.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, HashSet<Data>> result = cmsc.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, HashSet<Data>> result = cmqnRK.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, ArrayList<Data>> result = new_cpod4.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, ArrayList<Data>> result = new_cpod5.slide_process(incomingData, currentTime);
//          HashMap<OD_Query, ArrayList<Data>> result = new_cpod6.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, ArrayList<Data>> result = new_cpod7.slide_process(incomingData, currentTime);
//                HashMap<OD_Query, ArrayList<Data>> result = new_cpod8.slide_process(incomingData, currentTime);
                HashMap<OD_Query, ArrayList<Data>> result = new_cpod9.slide_process(incomingData, currentTime);
//                HashMap<OD_Query, ArrayList<Data>> result = new_cpod10.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, HashSet<Data>> result = new_cpod2.slide_process(incomingData, currentTime);
//            HashMap<OD_Query, ArrayList<Data>> result = mcsky.slide_process(incomingData, currentTime);
//           HashMap<OD_Query, ArrayList<Data>> result = sop.slide_process(incomingData, currentTime);

            double elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

            if (numberWindows > 1) {
                totalTime += elapsedTimeInSec;
            }

            for (OD_Query q : result.keySet()) {
                System.out.println("R = " + q.R + ", k = " + q.k + ", w = " + q.W + ", s = " + q.S + ", #outliers = " + result.get(q).size());
            }
            
            //compute avg outlier rate 
            int totalOutliers = 0;
            for(OD_Query q: result.keySet()){
                totalOutliers += result.get(q).size();
            }
            double cur_outlier_rate = totalOutliers*1.0/result.size()/Constants.W;
            overall_avg_outlier_rate = (overall_avg_outlier_rate * (numberWindows-1) + cur_outlier_rate)/(numberWindows);
            System.out.println("Overall avg outlier rate = "+ overall_avg_outlier_rate);
            System.out.println("Max memory = " + mesureThread.maxMemory * 1.0 / 1024 / 1024);
            mesureThread.averageTime = totalTime * 1.0 / (numberWindows - 1);
            mesureThread.writeResult();

            
        }

        mesureThread.averageTime = totalTime * 1.0 / (numberWindows - 1);
        mesureThread.writeResult();
        mesureThread.stop();
        mesureThread.interrupt();
        System.out.println("Num DCS = " + cmqnRK.numDCS / (numberWindows - 1));

        System.out.println("Time for processing expired slide = " + new_cpod.timeForProcessingExpiredSlide / new_cpod.count);
        System.out.println("Time for creating core = " + new_cpod.timeForCreatCore / new_cpod.count);
        System.out.println("Time for finding neighbors = " + new_cpod.timeForFindNeighbors / new_cpod.count);

    }

    public static double calculateAverage(List<Integer> marks) {
        Integer sum = 0;
        if (!marks.isEmpty()) {
            for (Integer mark : marks) {
                sum += mark;
            }
            return sum.doubleValue() / marks.size();
        }
        return sum;
    }

//    public static int getNewSlideSize()
    public static void readArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {

            //check if arg starts with --
            String arg = args[i];
            if (arg.indexOf("--") == 0) {
                switch (arg) {
                    case "--algorithm":
                        algorithm = args[i + 1];
                        break;
                    case "--R":
                        Constants.R = Double.valueOf(args[i + 1]);
                        break;
                    case "--W":
                        Constants.W = Integer.valueOf(args[i + 1]);
                        break;
                    case "--k":
                        Constants.k = Integer.valueOf(args[i + 1]);
                        Constants.minSizeOfCluster = Constants.k + 1;
                        break;
                    case "--datafile":
                        Constants.dataFile = args[i + 1];
                        break;
                    case "--output":
                        Constants.outputFile = args[i + 1];
                        break;
                    case "--numberWindow":
                        Constants.numberWindow = Integer.valueOf(args[i + 1]);
                        break;
                    case "--slide":
                        Constants.slide = Integer.valueOf(args[i + 1]);
                        break;
                    case "--resultFile":
                        Constants.resultFile = args[i + 1];
                        break;
                    case "--samplingTime":
                        Constants.samplingPeriod = Integer.valueOf(args[i + 1]);
                        break;
                    case "--matrixType":
                        Constants.matrixType = args[i + 1];
                        break;
                    case "--numCols":
                        Constants.numCols = Integer.valueOf(args[i + 1]);
                        break;

                }
            }
        }
    }

    public static void writeResult() {

        try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(Constants.resultFile, true)))) {
            for (Integer time : idOutliers) {
                out.println(time);
            }
        } catch (IOException e) {
        }

    }

}
