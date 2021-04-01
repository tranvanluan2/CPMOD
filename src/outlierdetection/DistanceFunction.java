/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import be.tarsos.lsh.Vector;
import com.sun.swing.internal.plaf.metal.resources.metal;
import java.util.Arrays;
import java.util.HashMap;
import mtree.tests.Data;
import mtree.tests.MTTest;
import outlierdetection.MCOD_UB_LB_New.MCO;

/**
 *
 * @author Luan
 */
public class DistanceFunction {

    // public static HashMap<String, Double> cache = new HashMap<>();
    public static double euclideanDistance(Data d1, Data d2) {
//        System.out.println("d1.values.length = "+d1.values.length);
//        System.out.println(Arrays.toString(d1.values));
        if(MTTest.count > 1){
            MTTest.numDCS +=1;
        }
        double sumSquare = 0;
        for (int i = 0; i < d1.values.length; i++) {
            sumSquare += Math.pow((d1.values[i] - d2.values[i]), 2);
        }
        double distance = Math.sqrt(sumSquare);
        return distance;
    }

    public static double euclideanDistance(Vector d1, Vector d2) {
//        System.out.println("d1.values.length = "+d1.values.length);
//        System.out.println(Arrays.toString(d1.values));
//        System.out.println("d.values.length = " + d1.values.length);
        double sumSquare = 0;
        for (int i = 0; i < d1.values.length; i++) {
            sumSquare += Math.pow((d1.values[i] - d2.values[i]), 2);
        }
        double distance = Math.sqrt(sumSquare);
        return distance;
    }

    public static double UB_distance1(MCOD_MESI_Overlap_Upper_Lower.MCData x,
            MCOD_MESI_Overlap_Upper_Lower.MCData y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] + y.std_list[i]) * (x.std_list[i] + y.std_list[i]));
        }
        return Math.sqrt(distance);

    }

    public static double UB_distance1(Vector x, Vector y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] + y.std_list[i]) * (x.std_list[i] + y.std_list[i]));
        }
        return Math.sqrt(distance);

    }

    public static double UB_distance1(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] + y.std_list[i]) * (x.std_list[i] + y.std_list[i]));
        }
        return Math.sqrt(distance);

    }

    static double UB_distance1(Upper_Lower_Mean_Distance3.C_Data x, Upper_Lower_Mean_Distance3.C_Data y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] + y.std_list[i]) * (x.std_list[i] + y.std_list[i]));
        }
        return Math.sqrt(distance);
    }

    public static double UB_distance2(Vector x, Vector y) {
        double distance = 0;
        for (int i = 0; i < x.values.length; i++) {
            distance += Math.abs(x.values[i] - y.values[i]);
        }
        return distance;
    }

    public static double UB_distance2(MCOD_MESI_Overlap_Upper_Lower.MCData x,
            MCOD_MESI_Overlap_Upper_Lower.MCData y) {
        double distance = 0;
        for (int i = 0; i < x.values.length; i++) {
            distance += Math.abs(x.values[i] - y.values[i]);
        }
        return distance;
    }

    public static double UB_distance2(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
        double distance = 0;
        for (int i = 0; i < x.values.length; i++) {
            distance += Math.abs(x.values[i] - y.values[i]);
        }
        return distance;
    }
    
    static double UB_distance2(Upper_Lower_Mean_Distance3.C_Data x, Upper_Lower_Mean_Distance3.C_Data y) {
        double distance = 0;
        for (int i = 0; i < x.values.length; i++) {
            distance += Math.abs(x.values[i] - y.values[i]);
        }
        return distance;
    }

    public static double UB_distance2(MCO x, MCO y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] + y.std_list[i]) * (x.std_list[i] + y.std_list[i]));
        }
        return Math.sqrt(distance);
    }

    public static double LB_distance2(Vector x, Vector y) {
        double distance = -1;
        int num_attribute = 3;
        for (int i = 0; i < num_attribute; i++) {
            distance = Math.max(distance, Math.abs(x.values[i] - y.values[i]));
        }
        return distance;
    }

    public static double LB_distance2(MCOD_MESI_Overlap_Upper_Lower.MCData x,
            MCOD_MESI_Overlap_Upper_Lower.MCData y) {
        double distance = -1;
        int num_attribute = 1;
        for (int i = 0; i < num_attribute; i++) {
            distance = Math.max(distance, Math.abs(x.values[i] - y.values[i]));
        }
        return distance;
    }

    public static double LB_distance2(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
        double distance = -1;
        int num_attribute = 1;
        for (int i = 0; i < num_attribute; i++) {
            distance = Math.max(distance, Math.abs(x.values[i] - y.values[i]));
        }
        return distance;
    }

    static double LB_distance2(Upper_Lower_Mean_Distance3.C_Data x, Upper_Lower_Mean_Distance3.C_Data y) {
        double distance = -1;
        for (int i: x.ordered_attributes.subList(0, 5)) {
            distance = Math.max(distance, Math.abs(x.values[i] - y.values[i]));
        }
        return distance;
    }

    public static double LB_distance1(MCOD_MESI_Overlap_Upper_Lower.MCData x,
            MCOD_MESI_Overlap_Upper_Lower.MCData y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] - y.std_list[i]) * (x.std_list[i] - y.std_list[i]));
        }

        return Math.sqrt(distance);
    }

    public static double LB_distance1(Vector x, Vector y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] - y.std_list[i]) * (x.std_list[i] - y.std_list[i]));
        }

        return Math.sqrt(distance);
    }

    public static double LB_distance1(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] - y.std_list[i]) * (x.std_list[i] - y.std_list[i]));
        }

        return Math.sqrt(distance);
    }

    static double LB_distance1(Upper_Lower_Mean_Distance3.C_Data x, Upper_Lower_Mean_Distance3.C_Data y) {
        double distance = 0;

        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] - y.std_list[i]) * (x.std_list[i] - y.std_list[i]));
        }

        return Math.sqrt(distance);

    }

    public static double LB_distance2(MCO x, MCO y) {
        double distance = 0;
        for (int i = 0; i < x.mean_list.length; i++) {
            distance += x.groups.get(i).length * ((x.mean_list[i] - y.mean_list[i]) * (x.mean_list[i] - y.mean_list[i])
                    + (x.std_list[i] - y.std_list[i]) * (x.std_list[i] - y.std_list[i]));
        }
        return Math.sqrt(distance);
    }

    public static int level_distance = 0;

//    public static double UB_distance(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
//        double distance = 0;
//        int d = x.dimensions();
//
//        distance += d * ((x.mean - y.mean) * (x.mean - y.mean) + (x.std + y.std) * (x.std + y.std));
//        double epsilon = 0.000000001;
//
////        System.out.println("x.std = "+x.std+", ");
//        for (int i = 0; i < level_distance; i++) {
//            distance -= x.std * y.std * (y.a[i] / (y.std + epsilon) + x.a[i] / (x.std + epsilon)) * (y.a[i] / (y.std + epsilon) + x.a[i] / (x.std + epsilon));
////            System.out.println("i = " + i +", distance = "+Math.sqrt(distance));
//
//        }
//        return Math.sqrt(distance);
//
//    }
//    public static double LB_distance(MCOD_MESI_Upper_Lower.MCData x, MCOD_MESI_Upper_Lower.MCData y) {
//        double distance = 0;
//        int d = x.dimensions();
//
//        distance += d * ((x.mean - y.mean) * (x.mean - y.mean) + (x.std - y.std) * (x.std - y.std));
//
//        double epsilon = 0.000000001;
//        for (int i = 0; i < level_distance; i++) {
//            distance += x.std * y.std * (y.a[i] / (y.std + epsilon) - x.a[i] / (x.std + epsilon)) * (y.a[i] / (y.std + epsilon) - x.a[i] / (x.std + epsilon));
//        }
//        return Math.sqrt(distance);
//    }

//    static double euclideanDistance(New_MQ_CPOD.C_Data d1, New_MQ_CPOD.CorePoint d2) {
//        double sumSquare = 0;
//        for (int i = 0; i < d1.values.length; i++) {
//            sumSquare += Math.pow((d1.values[i] - d2.orig_data.values[i]), 2);
//        }
//        double distance = Math.sqrt(sumSquare);
//        return distance;
//    }
//
//    static double euclideanDistance(New_MQ_CPOD.CorePoint d1, New_MQ_CPOD.CorePoint d2) {
//        double sumSquare = 0;
//        for (int i = 0; i < d1.orig_data.values.length; i++) {
//            sumSquare += Math.pow((d1.orig_data.values[i] - d2.orig_data.values[i]), 2);
//        }
//        double distance = Math.sqrt(sumSquare);
//        return distance;
//    }

    
}
