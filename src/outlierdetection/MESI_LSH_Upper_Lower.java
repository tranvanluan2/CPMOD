///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package outlierdetection;
//
//import be.tarsos.lsh.HashTable;
//import be.tarsos.lsh.Index;
//import be.tarsos.lsh.LSH;
//import be.tarsos.lsh.Vector;
//import be.tarsos.lsh.families.EuclidianHashFamily;
//import be.tarsos.lsh.families.HashFamily;
//import java.util.ArrayList;
//import java.util.Comparator;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.Hashtable;
//import java.util.List;
//import java.util.PriorityQueue;
//import mtree.tests.Data;
//import mtree.utils.Constants;
//import mtree.utils.Utils;
//
///**
// *
// * @author Luan Tran
// */
//public class MESI_LSH_Upper_Lower {
//
//    public static PriorityQueue<Vector> event_queue = new PriorityQueue(new VectorNeighborComparator());
//
//    public static int currentTime;
//
//    public static ArrayList<Double> p_distance_filtered_by_LB = new ArrayList<>();
//    public static ArrayList<Double> p_distance_filered_by_UB = new ArrayList<>();
//    public static long all_distance_computations = 0;
//    public static long filtered_b_LB = 0;
//
//    public static long all_distance_after_LB = 0;
//    public static long filtered_by_UB = 0;
//
//    public static long time_finding_neighbors = 0;
//    public static long count_data = 0;
//
//    HashFamily family;
//    public static HashMap<Integer, Index> all_lsh = new HashMap<>();
////    public static ArrayList<Integer> number_points_in_clusters = new ArrayList<>();
//
//    public MESI_LSH_Upper_Lower() {
//        if (Constants.dataFile.equals("covtype.data")) {
//            this.family = new EuclidianHashFamily(5000, 55);
//        }
//    }
//
//    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
//
//        ArrayList<Data> result = new ArrayList<>();
//        currentTime = _currentTime;
//        int expiredSlideIndex = (currentTime - Constants.W - 1) / slide;
////        System.out.println("Expire slide = "+ expiredSlideIndex);
////        System.out.println("LSH size "+ all_lsh.size());
////        System.out.println("Event Queue size "+ event_queue.size());
//        processExpiredData(expiredSlideIndex);
//
//        ArrayList<Vector> d_to_process = new ArrayList<>();
//
//        for (Data o : data) {
//            Vector d = new Vector(o);
//
//            d_to_process.add(d);
//
//            Index idx = all_lsh.get(d.getSlideIndex());
//            if (idx != null) {
//
//                idx.index(d);
//            } else {
//                idx = new Index(family, 3, 1);
//                idx.index(d);
//                all_lsh.put(d.getSlideIndex(), idx);
//            }
//
//        }
//
//        d_to_process.forEach((d) -> {
//            processData(d);
//        });
//
//        //find outlier     
//        for (Integer slideIndex : all_lsh.keySet()) {
//
//            List<Vector> ds = all_lsh.get(slideIndex).getAllData();
////            List<? super Vector> datas = ds;
//            for (Vector d : ds) {
//                int newestSlide = (currentTime - 1) / Constants.slide;
//                reProbe(d, newestSlide);
//                if (d.isOutlier()) {
//                    result.add(d);
//                }
//            }
//        }
//        return result;
//
//    }
//
//    public void reProbe(Vector d, int newestSlideIndex) {
//        int slideIndex = d.lastProbe + 1;
//        while (slideIndex <= newestSlideIndex) {
//            //find neighbor using LSH 
//
//            Index idx = all_lsh.get(slideIndex);
//
//            List<Vector> candidates = idx.query(d, Constants.k);
////                System.out.println("d.values.length = "+d.values.length);
////                System.out.println("Size of candidates = " + candidates.size());
//
//            System.out.println("************************************");
//            for (Vector d2 : candidates) {
////                    d2.values = d2_.values;
////                    checked.add(d2);
//                System.out.println("Distance = " + DistanceFunction.LB_distance2(d, d2));
//                if (d2.arrivalTime != d.arrivalTime && is_Neighbor2(d, d2)) {
//
//                    d.numSucceedingNeighbors += 1;
//                    d.countNeighbor += 1;
//
//                }
//                if (d.countNeighbor() >= Constants.k) {
//                    break;
//
//                }
//            }
//
//            if (!d.isOutlier()) {
//                d.lastProbe = slideIndex;
//                return;
//            }
//            d.lastProbe = slideIndex;
//            slideIndex++;
//        }
//    }
//
//    private void processExpiredData(int expiredSlide) {
//
//        all_lsh.remove(expiredSlide);
//
//        processEventQueue(expiredSlide);
//
//    }
//
//    private void processEventQueue(int expiredSlide) {
//        while (true) {
//            if (event_queue.isEmpty()) {
//                break;
//            }
//            Vector d = event_queue.peek();
//            if (d.getOldestSlide() <= expiredSlide || d.getSlideIndex() <= expiredSlide
//                    || d.numPrecedingNeighbor.isEmpty()) {
//                event_queue.poll();
//                if (d.numPrecedingNeighbor != null) {
//                    //d.updateEarliestNeighbor();
//                    for (int i = d.getOldestSlide(); i <= expiredSlide; i++) {
//                        if (d.numPrecedingNeighbor.get(i) != null) {
//
//                            d.countNeighbor -= d.numPrecedingNeighbor.get(i).size();
//
////                            d.numPrecedingNeighbor.get(i).clear();
//                            d.numPrecedingNeighbor.remove(i);
//                        }
//                    }
////                    d.updateEarliestNeighbor();
//                }
//                if (d.isUnsafeInlier() && !d.numPrecedingNeighbor.isEmpty() && d.getSlideIndex() > expiredSlide) {
//                    event_queue.add(d);
//                } else if (d.getSlideIndex() <= expiredSlide) {
//                    d.clean();
//                }
//            } else {
//                break;
//            }
//        }
//    }
//
//    private void processData(Vector d) {
//        long start = Utils.getCPUTime();
//        //find neighbor for d 
//        probe(d);
//        time_finding_neighbors += Utils.getCPUTime() - start;
//
//        if (d.isUnsafeInlier()) {
//            event_queue.add(d);
//        } else {
////            System.out.println("Not unsafe inlier, #succeeding neighbors = " + d.numSucceedingNeighbors);
////            System.out.println("Total neigbors = "+ d.countNeighbor());
//
//        }
//
//        count_data += 1;
//    }
//
//    public boolean is_Neighbor(Vector d, Vector d2) {
////        System.out.println("d.values.length = "+d.values.length);
////        System.out.println("LB distance = "+ DistanceFunction.LB_distance2(d, d2));
////        System.out.println("UB distance = "+ DistanceFunction.UB_distance2(d, d2));
////        System.out.println("Exact distance = "+ DistanceFunction.euclideanDistance(d, d2));
////        System.out.println("+++++++++++++++++++++++++++++");
//        if (DistanceFunction.LB_distance2(d, d2) <= Constants.R) {
//            return DistanceFunction.UB_distance2(d, d2) <= Constants.R
//                    || DistanceFunction.euclideanDistance(d, d2) <= Constants.R;
//        }
//        return false;
//
//    }
//
//    public boolean is_Neighbor2(Vector d, Vector d2) {
//
//        if (DistanceFunction.UB_distance2(d, d2) <= Constants.R) {
////            System.out.println("Good!");
//            return true;
//        } else {
//            return DistanceFunction.euclideanDistance(d, d2) <= Constants.R;
//        }
//
//    }
//
//    public void probe(Vector d) {
//
//        int slideIndex = d.getSlideIndex();
//
//        while (d.isOutlier() && (slideIndex > d.getSlideIndex() - Constants.W / Constants.slide) && slideIndex >= 0) {
//
////            System.out.println("d.values.length = "+d.values.length);
//            ArrayList<Vector> neighbors = new ArrayList<>();
////            HashSet<Vector> checked = new HashSet<>(Constants.slide);
//            //find neighbor using LSH 
////            if (all_lsh.containsKey(slideIndex)) {
////            System.out.println("Slide index =" + slideIndex);
////            System.out.println("Key lsh = ");
////            for(Integer k: all_lsh.keySet()){
////                System.out.print(k+",");
////            }
//            Index idx = all_lsh.get(slideIndex);
//
////            
////            for(HashTable h: idx.hashTable){
////                for(Integer k: h.hashTable.keySet()){
////                    System.out.println("Key = "+ k);
////                }
////                System.out.println("+++++++++++++++++++++++++++++++++++++++");
////            }
//            List<Vector> candidates = idx.query(d, Constants.k);
////                System.out.println("d.values.length = "+d.values.length);
////                System.out.println("Size of candidates = " + candidates.size());
//            int count = 0;
//            for (Vector d2 : candidates) {
////                    d2.values = d2_.values;
////                    checked.add(d2);
//                count += 1;
//                if (d2.arrivalTime != d.arrivalTime && is_Neighbor(d, d2)) {
//                    if (d.getSlideIndex() <= slideIndex) {
//                        d.numSucceedingNeighbors += 1;
//                        d.countNeighbor += 1;
//                    } else {
//                        neighbors.add(d2);
//                    }
//
//                }
//                if (d.countNeighbor() + neighbors.size() >= Constants.k) {
//                    break;
//
//                }
//            }
//
//            if (slideIndex < d.getSlideIndex()) {
//                d.numPrecedingNeighbor.put(slideIndex, neighbors);
//                d.countNeighbor += neighbors.size();
//            }
//
//            slideIndex--;
//        }
//
//    }
//
//    private static class VectorNeighborComparator implements Comparator<Vector> {
//
//        @Override
//        public int compare(Vector t, Vector t1) {
//            if (t.getOldestSlide() > t1.getOldestSlide()) {
//                return 1;
//            } else {
//                return -1;
//            }
//        }
//    }
//
//}
