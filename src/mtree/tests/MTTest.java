package mtree.tests;

import be.tarsos.lsh.Vector;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import outlierdetection.AbstractC;
import outlierdetection.ApproxStorm;
import outlierdetection.Direct_Update_Event;
import outlierdetection.ExactStorm;
import outlierdetection.Lazy_Update_Event;
import outlierdetection.MESI;
import outlierdetection.MicroCluster;
import mtree.utils.Constants;
import mtree.utils.Utils;
import outlierdetection.ApproxCPOD2;
import outlierdetection.DataLUEObject;
import outlierdetection.IMCOD;
import outlierdetection.MCOD_MESI;
import outlierdetection.MCOD_MESI_Overlap_Upper_Lower;
import outlierdetection.MCOD_MESI_Safe;
import outlierdetection.MCOD_MESI_Upper_Lower;
import outlierdetection.MCOD_Safe_Version;
import outlierdetection.MCOD_Safe_Wait;
import outlierdetection.MCOD_UB_LB_New;
import outlierdetection.MCOD_UP;
import outlierdetection.MESIWithHash;
//import outlierdetection.MESI_LSH_Upper_Lower;
import outlierdetection.MesiAndCluster;
import outlierdetection.MicroCluster_New;
import outlierdetection.MicroCluster_NewVersion;
import outlierdetection.MicroCluster_NewVersion.MCData;
import outlierdetection.NewCorePoint;
import outlierdetection.UpperLower;
import outlierdetection.Upper_Lower_Mean_Distance;
import outlierdetection.Upper_Lower_Mean_Distance3;
import outlierdetection.Upper_Lower_Mean_Distance4;
import outlierdetection.Upper_Lower_Mean_Distance_2;
import outlierdetection.Upper_Lower_Mean_Distance_5;

public class MTTest {

    public static int currentTime = 0;

    public static boolean stop = false;

    public static HashSet<Integer> idOutliers = new HashSet<>();

    public static String algorithm;
    public static int numberWindows = 0;
    public static double average_number_cluster = 0;

    public static double[] max_values;
    public static double[] min_values;

    public static double start;

    public static double timeForIndexing = 0;
    public static double timeForNeighborSearch = 0;
    public static long numDCS = 0;
    public static int count = 0;
//    public static long numPointsNeedNS = 0;
//    public static long numDSForNS = 0;
//    public static long numDSForIndex = 0;
    public static void main(String[] args) {

        readArguments(args);
        for (String arg : args) {
            System.out.print(arg + " ");
        }
        System.out.println("");

        MesureMemoryThread mesureThread = new MesureMemoryThread();
        mesureThread.start();
//         Stream s = Stream.getInstance("ForestCover");
        Stream s = Stream.getInstance("");
//         Stream s = Stream.getInstance("randomData");
//        Stream s = Stream.getInstance("randomData1");
        // Stream s = Stream.getInstance(null);
        // Stream s = Stream.getInstance("tagData");
//        Stream s = Stream.getInstance("Trade");

        ExactStorm estorm = new ExactStorm();
        ApproxStorm apStorm = new ApproxStorm(0.1);
        AbstractC abstractC = new AbstractC();
        Lazy_Update_Event lue = new Lazy_Update_Event();
        Direct_Update_Event due = new Direct_Update_Event();
        MicroCluster micro = new MicroCluster();
        MicroCluster_New mcnew = new MicroCluster_New();
        MesiAndCluster mac = new MesiAndCluster();
        MESI mesi = new MESI();
        MESIWithHash mesiWithHash = new MESIWithHash();
        IMCOD imcod = new IMCOD();
        MicroCluster_NewVersion mcod_new = new MicroCluster_NewVersion();
        MCOD_Safe_Version mcod_safe = new MCOD_Safe_Version();
        MCOD_MESI_Safe mcod_mesi_safe = new MCOD_MESI_Safe();
        MCOD_Safe_Wait mcod_safe_wait = new MCOD_Safe_Wait();
        MCOD_MESI mcod_mesi = new MCOD_MESI();
        UpperLower upl = new UpperLower();
        MCOD_UP mcod_up = new MCOD_UP();
        MCOD_MESI_Upper_Lower mcod_up_lo = new MCOD_MESI_Upper_Lower();
        MCOD_MESI_Overlap_Upper_Lower mcod_mesi_overlap_ub_lb = new MCOD_MESI_Overlap_Upper_Lower();
//        MESI_LSH_Upper_Lower mesi_ul = new MESI_LSH_Upper_Lower();
        Upper_Lower_Mean_Distance upper_lower_distance = new Upper_Lower_Mean_Distance();
        Upper_Lower_Mean_Distance_2 upper_lower_distance_2 = new Upper_Lower_Mean_Distance_2();
        Upper_Lower_Mean_Distance3 upper_lower_distance_3 = new Upper_Lower_Mean_Distance3();
        Upper_Lower_Mean_Distance4 upper_lower_distance_4 = new Upper_Lower_Mean_Distance4();
        Upper_Lower_Mean_Distance_5 upper_lower_distance_5 = new Upper_Lower_Mean_Distance_5();
        ApproxCPOD2 approx_cpod = new ApproxCPOD2();
        MCOD_UB_LB_New mcod_ub_lb = new MCOD_UB_LB_New();
        NewCorePoint ncp = new NewCorePoint();
//        Approx ncp = new NewCorePoint();

        double avg_core = 0;
        
        double totalTime = 0;
        File file = new File("/home/luan/Windows_Data2/MTree_New/detected_outlier_idxs.txt"); 
        if(file.delete()) 
        { 
            System.out.println("File deleted successfully"); 
        } 
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
                
                
//                System.out.println("Last idx time = " + (incomingData.get(incomingData.size()-1).arrivalTime-1));
                //init max, min values
                Constants.dimensions = incomingData.get(0).values.length;
                max_values = new double[Constants.dimensions];
                min_values = new double[Constants.dimensions];
                //update max_values, min_values
                for (int i = 0; i < Constants.dimensions; i++) {
                    max_values[i] = incomingData.get(0).values[i];
                    min_values[i] = incomingData.get(0).values[i];
                    for (Data d : incomingData) {
                        if (d.values[i] > max_values[i]) {
                            max_values[i] = d.values[i];
                        }
                        if (d.values[i] < min_values[i]) {
                            min_values[i] = d.values[i];
                        }
                    }
                }

            }

//            //normalize incoming data
//            for (int i = 0; i < Constants.dimensions; i++) {
//
//                for (Data d : incomingData) {
//                    d.values[i] = (d.values[i] - min_values[i]) / (max_values[i] - min_values[i]);
//                }
//            }
            start = Utils.getCPUTime(); // requires java 1.5

            /**
             * do algorithm
             */
            switch (algorithm) {
//                case "mesi_upper_lower":
//                    ArrayList<Data> outliers = mesi_ul.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
//                    System.out.println("Num outliers = " + outliers.size());
//                    double elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//
//                    if (numberWindows > 1) {
//                        totalTime += elapsedTimeInSec;
//                    }
//                    
////                    outliers.stream().forEach((outlier) -> {
////                        idOutliers.add(outlier.arrivalTime);
////                    });
//                    
//                    break;

                case "mcod_mesi_upper_lower":
                    ArrayList<Data> outliers = mcod_up_lo.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    double elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                    
                case "approx_cpod":
                    outliers = approx_cpod.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);

                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    
                    System.out.println("Num outliers = " + outliers.size());
                    break;
                case "upper_lower_mean_distance5":
                    outliers = upper_lower_distance_5.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);

                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    System.out.println("Num outliers = " + outliers.size());
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    int max_id_outlier = 0;
//                    int min_id_outlier = Integer.MAX_VALUE;
                    //write the outlier arrival time to file
//                    try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(Constants.outputFile, true)))) {
//                        for (Data d : outliers) {
//                            out.println(d.arrivalTime-1);
//                            if(d.arrivalTime-1> max_id_outlier){
//                                max_id_outlier = d.arrivalTime-1;
//                            }
//                            if(d.arrivalTime-1 < min_id_outlier){
//                                min_id_outlier = d.arrivalTime-1;
//                            }
//                        }
//                        System.out.println("Max id outlier = "+ max_id_outlier);
//                        System.out.println("Min id outlier = "+ min_id_outlier);
//                    } catch (IOException e) {
//                    }
                    System.out.println("Elapsed Time in Sec = "+ elapsedTimeInSec);

                    if (count > 0) {
                        System.out.println("Average #DCS = "+ numDCS/count);
                    }
                    avg_core =(avg_core * (numberWindows-1)+Upper_Lower_Mean_Distance_5.all_distinct_cores.size() )/numberWindows;
                    System.out.println("Num distinct cores = "+avg_core );
                    System.out.println("ALL DISTINCT CORES = "+ Upper_Lower_Mean_Distance_5.all_distinct_cores.size());
//                    System.out.println("Num points need NS = "+ Upper_Lower_Mean_Distance_5.numPointNeedNS*1.0/numberWindows);
//                    System.out.println("Avg num DCS =" + Upper_Lower_Mean_Distance_5.numDCS*1.0/Upper_Lower_Mean_Distance_5.numPointNeedNS);
//                    System.out.println("Avg total for each w =" + 
//                            (Upper_Lower_Mean_Distance_5.numDCS +Upper_Lower_Mean_Distance_5.numDCSForIndexing) *1.0/numberWindows);
                    break;

                case "upper_lower_mean_distance4":
                    outliers = upper_lower_distance_4.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                case "upper_lower_mean_distance3":
                    outliers = upper_lower_distance_3.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;

                case "new_core_point":
                    outliers = ncp.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                case "mcod_mesi_overlap_ub_lb":
                    outliers = mcod_mesi_overlap_ub_lb.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;

                case "mcod_ub_lb":
                    outliers = mcod_ub_lb.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;

                case "upper_lower_mean_distance":
                    outliers = upper_lower_distance.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
//                    if(numberWindows < 100){
//                        upper_lower_distance.slide_range =1;
//                    }
//                    else{
//                        upper_lower_distance.slide_range = 1;
//                    }

//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
//                    System.out.println("Time for processing expired slide  = "+
//                            Upper_Lower_Mean_Distance.timeForProcessingExpiredSlide/numberWindows);
//                    System.out.println("Time for new data = "+ Upper_Lower_Mean_Distance.timeForProcessingNewData/numberWindows);
//                    System.out.println("Time for reprobing data = "+ Upper_Lower_Mean_Distance.timeForProcessingReprobing/numberWindows);
//                    System.out.println("Total Check neighbors = "+Upper_Lower_Mean_Distance.countCheckNeighbor/numberWindows);
//                    System.out.println("Total Filted by LB neighbors = "+Upper_Lower_Mean_Distance.countFiltedByLB/numberWindows);
                    break;

                case "upper_lower_mean_distance_2":
                    outliers = upper_lower_distance_2.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    break;

                case "mcod_up":

                    outliers = mcod_up.detectOutlier(incomingData);
                    System.out.println("Num outliers = " + outliers.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }

                    outliers.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    break;
                case "upper_lower":
                    outliers = upl.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    break;
                case "exactStorm":
                    outliers = estorm.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
//                    outliers.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });

                    break;
                case "compareApproxExact":
//                    ArrayList<MCData> exactResult = mcod_new.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
//                    ArrayList<Data> approxResult = apStorm.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
//                    if (numberWindows > 1) {
//                        computePrecesionRecall(exactResult, approxResult);
//                        System.out.println("Precision = " + precision);
//                        System.out.println("Recall = " + recall);
//                    }

                case "mesiAndCluster":
                    ArrayList<Data> outliers100 = mac.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    System.out.println("Num outlier = " + outliers100.size());
                    totalTime += elapsedTimeInSec;
//                    outliers100.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });

                    break;
                case "imcod":
                    ArrayList<Data> outliers200 = imcod.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    System.out.println("Num outliers = " + outliers200.size());
//                    outliers200.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                case "mcod_new":
                    ArrayList<MCData> outliers300 = mcod_new.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    System.out.println("num outliers = " + outliers300.size());
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    outliers300.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    average_number_cluster += MicroCluster_NewVersion.microClusters.size();
                    break;
                case "mcod_mesi":
                    ArrayList<Data> outliers123 = mcod_mesi.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    System.out.println("Num outliers = " + outliers123.size());
                    System.out.println("Percentage added to cluster = " + MCOD_MESI.p_add_to_cluster);
                    if (count > 0) {
                        System.out.println("Average #DCS = "+ numDCS/count);
                    }
                     System.out.println("Elapsed Time in Sec = "+ elapsedTimeInSec);
//                    System.out.println("Num points need NS "+ MTTest.numPointsNeedNS/MTTest.numberWindows);
//                    System.out.println("AVG DCS per points "+ MTTest.numDSForNS*1.0/MTTest.numPointsNeedNS);
//                    System.out.println("Total DCS each slide =" + (MTTest.numDSForIndex + MTTest.numDSForNS)*1.0/MTTest.numberWindows);
//                    outliers123.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
//                    average_number_cluster += MCOD_MESI.microClusters.size();
                    break;

                case "mcod_safe":
                    ArrayList<MCOD_Safe_Version.MCData> outliers400 = mcod_safe.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    outliers400.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    average_number_cluster += MCOD_Safe_Version.safeClusters.size();
                    break;
                case "mcod_safe_wait":
                    ArrayList<Data> outliers_safe_wait = mcod_safe_wait.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    outliers_safe_wait.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    break;
                case "mcod_mesi_safe":
                    ArrayList<MCOD_MESI_Safe.MCData> outliers500 = mcod_mesi_safe.detectOutlier(incomingData, currentTime, Constants.W, Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
                    System.out.println("Num outliers = " + outliers500.size());
                    outliers500.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    break;
                case "approximateStorm":
                    ArrayList<Data> outliers2 = apStorm.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
//                    outliers2.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                case "abstractC":
                    ArrayList<Data> outliers3 = abstractC.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
                    outliers3.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);

                    });
                    break;
                case "lue":
                    HashSet<DataLUEObject> outliers4 = lue.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
//                    outliers4.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    break;
                case "due":
                    HashSet<DataLUEObject> outliers5 = due.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
                    outliers5.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    break;
                case "microCluster":
                    ArrayList<Data> outliers6 = micro.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    outliers6.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);

                    });

                    break;
                case "microCluster_new":
                    ArrayList<Data> outliers9 = mcnew.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                    if (numberWindows > 1) {
                        totalTime += elapsedTimeInSec;
                    }
                    outliers9.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);

                    });

                    System.out.println("Num outliers = " + outliers9.size());
                    if (count > 0) {
                        System.out.println("Average #DCS = "+ numDCS/count);
                    }
                    
                     System.out.println("Elapsed Time in Sec = "+ elapsedTimeInSec);

//                    ArrayList<Data> outliers10 = estorm.detectOutlier(incomingData, currentTime, Constants.W,
//                            Constants.slide);
//                    
//                    System.out.println("--------------------------------------------------");
//                    System.out.println("Not in exact storm");
//                    for(Data d: outliers9){
//                        if(!outliers10.contains(d))
//                            System.out.println(d.arrivalTime);
//                    }
//                    System.out.println("---------------------------------------------------");
                    break;
                case "mesi":
                    ArrayList<Data> outliers7 = mesi.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
//                    outliers7.stream().forEach((outlier) -> {
//                        idOutliers.add(outlier.arrivalTime);
//                    });
                    System.out.println("Num outlier = " + outliers7.size());
                    break;
                case "mesiWithHash":
                    HashSet<Vector> outliers8 = mesiWithHash.detectOutlier(incomingData, currentTime, Constants.W,
                            Constants.slide);
                    elapsedTimeInSec = (Utils.getCPUTime() - start) * 1.0 / 1000000000;

                    totalTime += elapsedTimeInSec;
                    outliers8.stream().forEach((outlier) -> {
                        idOutliers.add(outlier.arrivalTime);
                    });
                    break;
                    
                 

            }
            
            if (numberWindows == 1) {
                totalTime = 0;
                MesureMemoryThread.timeForIndexing = 0;
                MesureMemoryThread.timeForNewSlide = 0;
                MesureMemoryThread.timeForExpireSlide = 0;
                MesureMemoryThread.timeForQuerying = 0;

//                MicroCluster_New.timeForAddToCluster = 0;
//                MicroCluster_New.timeForAddToPD = 0;
//                //MicroCluster_New. = 0;
//                MicroCluster_New.timeForFindCluster = 0;
//                MicroCluster_New.timeForFormNewCluster = 0;
//                MicroCluster_New.timeForRemovePointFromCluster = 0;
//                MicroCluster_New.timeForRemovePointFromPD = 0;
//                MicroCluster_New.timeForUpdateAffectedPointInCluster = 0;
//                MicroCluster_New.timeForUpdateAffectedPointInEventQueue = 0;
            }
//            System.out.println("#window: " + numberWindows);
//            System.out.println("Total #outliers: " + idOutliers.size());
//            System.out.println("Average Time: " + totalTime * 1.0 / numberWindows);
//            System.out.println("Peak memory: " + MesureMemoryThread.maxMemory * 1.0 / 1024 / 1024);
//            System.out.println("Time index, remove data from structure: " + MesureMemoryThread.timeForIndexing * 1.0 / 1000000000 / numberWindows);
//            System.out.println("Time for querying: " + MesureMemoryThread.timeForQuerying * 1.0 / 1000000000 / numberWindows);
//            System.out.println("Time for new slide: " + MesureMemoryThread.timeForNewSlide * 1.0 / 1000000000 / numberWindows);
//            System.out.println("Time for expired slide: " + MesureMemoryThread.timeForExpireSlide * 1.0 / 1000000000 / numberWindows);
//            System.out.println("------------------------------------");
//
//            if (algorithm.equals("exactStorm")) {
//
//                System.out.println("Avg neighbor list length = " + ExactStorm.avgAllWindowNeighbor / numberWindows);
//            } else if (algorithm.equals("mesi")) {
//
//                System.out.println("Avg trigger list = " + MESI.avgAllWindowTriggerList / numberWindows);
//                System.out.println("Avg neighbor list = " + MESI.avgAllWindowNeighborList / numberWindows);
//            } else if (algorithm.equals("microCluster")) {
//
//                System.out.println("Number clusters = " + MicroCluster.numberCluster / numberWindows);
//                System.out.println("Max  Number points in event queue = " + MicroCluster.numberPointsInEventQueue);
//
//                System.out.println("Avg number points in clusters= " + MicroCluster.numberPointsInClustersAllWindows / numberWindows);
//                System.out.println("Avg Rmc size = " + MicroCluster.avgPointsInRmcAllWindows / numberWindows);
//                System.out.println("Avg Length exps= " + MicroCluster.avgLengthExpsAllWindows / numberWindows);
//            } else if (algorithm.equals("due")) {
////            Direct_Update_Event.numberPointsInEventQueue = Direct_Update_Event.numberPointsInEventQueue /numberWindows;
//                Direct_Update_Event.avgAllWindowNumberPoints = Direct_Update_Event.numberPointsInEventQueue;
//                System.out.println("max #points in event queue = " + Direct_Update_Event.avgAllWindowNumberPoints);
//            }
//            if (algorithm.equals("microCluster_new")) {
//                System.out.println("avg points in clusters = "+MicroCluster_New.avgNumPointsInClusters *1.0/numberWindows);
//                System.out.println("Avg points in event queue = "+ MicroCluster_New.avgNumPointsInEventQueue*1.0/numberWindows);
//                System.out.println("avg neighbor list length = "+ MicroCluster_New.avgNeighborListLength*1.0/numberWindows);
//                System.out.println("Time for forming new cluster = "+ MicroCluster_New.timeForFormNewCluster*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for adding to cluster= "+ MicroCluster_New.timeForAddToCluster*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for adding to pd= "+ MicroCluster_New.timeForAddToPD*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for remove from cluster= "+ MicroCluster_New.timeForRemovePointFromCluster*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for remove from pd= "+ MicroCluster_New.timeForRemovePointFromPD*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for updating affected points in cluster= "+ MicroCluster_New.timeForUpdateAffectedPointInCluster*1.0/numberWindows/ 1000000000 );
//                System.out.println("Time for updating affected points in pd= "+ MicroCluster_New.timeForUpdateAffectedPointInEventQueue*1.0/numberWindows/ 1000000000 );
//
//            }

        }

//       
////        Constants.numberWindow--;
//        ExactStorm.avgAllWindowNeighbor = ExactStorm.avgAllWindowNeighbor / numberWindows;
//        MESI.avgAllWindowTriggerList = MESI.avgAllWindowTriggerList / numberWindows;
//        MicroCluster.numberCluster = MicroCluster.numberCluster / numberWindows;
//        MicroCluster.avgPointsInRmcAllWindows = MicroCluster.avgPointsInRmcAllWindows / numberWindows;
//        MicroCluster.avgLengthExpsAllWindows = MicroCluster.avgLengthExpsAllWindows / numberWindows;
//        MicroCluster.numberPointsInClustersAllWindows = MicroCluster.numberPointsInClustersAllWindows / numberWindows;
//        MicroCluster_New.avgNumPointsInClusters = MicroCluster_New.avgNumPointsInClusters/numberWindows;
        mesureThread.averageTime = totalTime * 1.0 / (numberWindows - 1);
        mesureThread.writeResult();
        mesureThread.stop();
        mesureThread.interrupt();

        System.out.println("Avg time for indexing = "+ timeForIndexing*1.0/(numberWindows-1));
        System.out.println("Avg time for neighbor search = "+ timeForNeighborSearch*1.0/(numberWindows-1));
        
        System.out.println("Write core points info");
//        System.out.println("Average numer of clusters = "+ average_number_cluster*1.0/numberWindows);
//        System.out.println("Average number of point in clusters = "+ MicroCluster_NewVersion.avgPointsInClusters);
        /**
         * Write result to file
         */
//        System.out.println("Total number of outliers = "+ idOutliers.size());
//        if (!"".equals(Constants.resultFile)) {
//            writeResult();
//        }
        //print statistics
//        double total_p_added_to_cluster = 0;
//        for(Double p: MCOD_MESI_Upper_Lower.p_add_to_cluster){
//            total_p_added_to_cluster += p;
//        }
//        System.out.println("Percentage added to cluster = "+ total_p_added_to_cluster/MCOD_MESI_Upper_Lower.p_add_to_cluster.size());
//        System.out.println("Percentage filtered by LB = "+ MCOD_MESI_Upper_Lower.filtered_b_LB*1.0/MCOD_MESI_Upper_Lower.all_distance_computations);
//        System.out.println("Percentage filtered by UB = "+ MCOD_MESI_Upper_Lower.filtered_by_UB*1.0/MCOD_MESI_Upper_Lower.all_distance_after_LB);
////      
//        System.out.println("Time for finding clusters = "+ MCOD_MESI_Upper_Lower.time_adding_to_cluster*1.0/MCOD_MESI_Upper_Lower.count_data/1000000000 );
//        System.out.println("Time for finding neighbors  = "+ MCOD_MESI_Upper_Lower.time_finding_neighbors*1.0/MCOD_MESI_Upper_Lower.count_data/1000000000 );
//        
//        System.out.println("Time for finding clusters = "+ mcod_new.time_adding_to_cluster*1.0/MCOD_MESI_Upper_Lower.count_data/1000000000 );
//        System.out.println("Time for finding neighbors  = "+ MCOD_MESI_Upper_Lower.time_finding_neighbors*1.0/MCOD_MESI_Upper_Lower.count_data/1000000000 );
    }

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

    public static double precision = 0.0;
    public static double recall = 0.0;

    public static int countPrecisionRecalll = 0;

    private static void computePrecesionRecall(ArrayList<Data> exactResult, ArrayList<Data> approxResult) {
//        HashSet<Integer> exact = new HashSet<>();
//        exactResult.stream().forEach((d) -> {
//            exact.add(d.arrivalTime);
//        });
//
//        HashSet<Integer> approx = new HashSet<>();
//        approxResult.stream().forEach((d) -> {
//            approx.add(d.arrivalTime);
//        });
        //compute precision
        countPrecisionRecalll++;
        int correctDetect = 0;
        for (Data d : approxResult) {
            correctDetect = exactResult.stream().filter((d2) -> (d2.arrivalTime == d.arrivalTime)).map((_item) -> 1).reduce(correctDetect, Integer::sum);
        }
        if (!exactResult.isEmpty()) {
            double newPrecision = correctDetect * 1.0 / exactResult.size();
            //update precision
            if (newPrecision != Double.NaN) {
                precision = (precision * (countPrecisionRecalll - 1) + newPrecision) / countPrecisionRecalll;
            }
        } else if (approxResult.isEmpty()) {
            double newPrecision = 1;
            //update precision
            if (newPrecision != Double.NaN) {
                precision = (precision * (countPrecisionRecalll - 1) + newPrecision) / countPrecisionRecalll;
            }
        } else {
            double newPrecision = 1;
            //update precision
            if (newPrecision != Double.NaN) {
                precision = (precision * (countPrecisionRecalll - 1) + newPrecision) / countPrecisionRecalll;
            }
        }
        if (!approxResult.isEmpty()) {
            //compute recall 
            double newRecall = correctDetect * 1.0 / approxResult.size();
            //update recall 
            if (newRecall != Double.NaN) {
                recall = (recall * (countPrecisionRecalll - 1) + newRecall) / countPrecisionRecalll;
            }
        } else {
            //compute recall 
            double newRecall = 1;
            //update recall 
            if (newRecall != Double.NaN) {
                recall = (recall * (countPrecisionRecalll - 1) + newRecall) / countPrecisionRecalll;
            }
        }
    }
}
