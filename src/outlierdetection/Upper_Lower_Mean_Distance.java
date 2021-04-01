/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import be.tarsos.lsh.Vector;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author Luan Tran
 */
public class Upper_Lower_Mean_Distance {

//    public static PriorityQueue<Vector> event_queue = new PriorityQueue(new VectorNeighborComparator());
    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<Vector>> all_slides = new HashMap<>();
//    public static ArrayList<Integer> number_points_in_clusters = new ArrayList<>();

    public static HashMap<Integer, ArrayList<Vector>> all_core_points = new HashMap<>();

//    public static HashMap<Integer, HashSet<Vector>> trigger_list = new HashMap<>();
    public static HashMap<Integer, HashSet<Vector>> outlierList = new HashMap<>();

    public static int countCheckNeighbor = 0;
    public static int countFiltedByLB = 0;
//    public static int countFiltedByUB = 0;

    public static double timeForProcessingExpiredSlide = 0;
    public static double timeForProcessingNewData = 0;
    public static double timeForProcessingReprobing = 0;

    public static HashMap<Integer, HashSet<Vector>> corePointMapTrigger = new HashMap<>();
    public static HashMap<Integer, HashSet<Vector>> neighborCountTrigger = new HashMap<>();

    public static int numFoundCore = 0;

    double avg_candidate_size = 0;
    int countCore = 0;
    double avg_closeNeighbor = 0;

    public static double avg_depth = 0;
    public static int countProb = 0;

    public Upper_Lower_Mean_Distance() {

    }

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
        numFoundCore = 0;
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        currentTime = _currentTime;
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W/slide;
        System.out.println("Expire slide index = " + expiredSlideIndex);
        
        
        long start = Utils.getCPUTime();
        processExpiredData(expiredSlideIndex);
        timeForProcessingExpiredSlide += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        System.out.println("Time for processing expired slide = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        
        ArrayList<Vector> d_to_process = new ArrayList<>(data.size());
//        HashSet<Integer> slide_to_process = new HashSet<>();

        int[] slide_to_process;
        if (data.size() == W) {
            slide_to_process = new int[W / slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
        } else {
            slide_to_process = new int[]{(currentTime - 1) / slide};
        }
//        for(Data o: data)
        for (int i = 0; i < data.size(); i++) {
            Data o = data.get(i);
            Vector d = new Vector(o);
//            Random r = new Random();
//            if (r.nextInt(1000) < 500) {
//
//                //check if d is not neighbor of any existing core point
//                boolean isNeighbor = false;
//                ArrayList<Vector> all_existing_cores = all_core_points.get(d.getSlideIndex());
//                if (all_existing_cores != null) {
//                    for (Vector v : all_existing_cores) {
//                        if (check_distance_neighbor(d, v) <= Constants.R / 2) {
//                            isNeighbor = true;
//                            break;
//                        }
//                    }
//                }
//                if (!isNeighbor) {
//                    d.isCore = true;
//                    if (all_core_points.containsKey(d.getSlideIndex())) {
//                        all_core_points.get(d.getSlideIndex()).add(d);
//                    } else {
//                        ArrayList<Vector> c = new ArrayList<>();
//                        c.add(d);
//                        all_core_points.put(d.getSlideIndex(), c);
//                    }
//                }
//            }

            d_to_process.add(d);

            ArrayList<Vector> idx = all_slides.get(d.getSlideIndex());
            if (idx != null) {

                idx.add(d);
            } else {
                idx = new ArrayList<>(Constants.slide);
                idx.add(d);
                all_slides.put(d.getSlideIndex(), idx);
            }
        }

        start = Utils.getCPUTime();
        for (Integer sIdx : slide_to_process) {

            ArrayList<Vector> corePoints = selectCore(sIdx);
            all_core_points.put(sIdx, corePoints);
            System.out.println("Num new core = " + corePoints.size());
//            for (Vector c : corePoints) {
//                System.out.println("Num close neighbors = " + c.closeNeighbors.size());
//            }
        }
        System.out.println("Time for selecting new core = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        //process core points 
        if (!all_core_points.isEmpty()) {
            for (int s : slide_to_process) {
                if (all_core_points.get(s) != null) {
                    for (Vector c : all_core_points.get(s)) {
//                        prob_for_core(c);
                        avg_candidate_size += c.doubleRCandidates.size();
                        avg_closeNeighbor += c.closeNeighbors.size();
                        countCore += 1;
                    }
                }
            }
        }
        System.out.println("AVG Candidate double size = " + avg_candidate_size / countCore);
        System.out.println("AVG close neighbor size = " + avg_closeNeighbor / countCore);
        start = Utils.getCPUTime();
        int newestSlide = (currentTime - 1) / Constants.slide;
        for (int i = 0; i < d_to_process.size(); i++) {
//        for(Vector d: d_to_process){
            Vector d = d_to_process.get(i);
            if (d.isOutlier()) {
                probe(d, newestSlide, false);
            }
        }
        System.out.println("Avg depth = " + avg_depth);
        System.out.println("Time for probing new data = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//        System.out.println("Num Found Core = " + numFoundCore);
        timeForProcessingNewData += (Utils.getCPUTime() - start) * 1.0 / 1000000000;

        start = Utils.getCPUTime();

//        for (int slideIndex : all_slides.keySet()) {
//
//            List<Vector> ds = all_slides.get(slideIndex);
//            for (int i = 0; i < ds.size(); i++) {
////            for(Vector d: ds){
//                Vector d = ds.get(i);
//                if (d.isOutlier() && d.lastProb < newestSlide) {
//
//                    reProbe(d, newestSlide);
//                }
//                if (d.isOutlier()) {
//                    result.add(d);
//                }
//            }
//
//        }
        for (int slideIndex : outlierList.keySet()) {
            for (Vector d : outlierList.get(slideIndex)) {
                if (d.isOutlier()) {
                    if (d.lastProbRight == -1) {
                        reProbe(d, newestSlide);
                        if (d.isOutlier()) {
                            probe(d, d.getSlideIndex(), false);
                        }
                    } else if (d.lastProbRight < newestSlide) {
                        reProbe(d, newestSlide);
                    }

                }
                if (d.isOutlier()) {
                    result.add(d);
                }
            }
        }
        timeForProcessingReprobing += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        System.out.println("Time for reprobing = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        return result;

    }

    public void reProbe(Vector d, int newestSlideIndex) {
//        d.neighborCount -= d.numNeighborPassive;
//        d.numNeighborPassive = 0;
//  find safe core point in succeeding slides
//        int count = 0;
//        for (int slideIndex = d.lastProb + 1; slideIndex <= newestSlideIndex; slideIndex++) {
//            if (all_core_points.containsKey(slideIndex)) {
//                for (Vector c : all_core_points.get(slideIndex)) {
//                    if (check_distance_R2neighbor(d, c) <= Constants.R / 2) {
//                        d.neighborCount += c.closeNeighbors.size();
//                        count += c.closeNeighbors.size();
//                        if (d.neighborCount >= Constants.k) {
////                            System.out.println("Found!");
//                            return;
//                        }
////                    System.out.println("Found!!!!!!!!!!!!!!!");
////                    return;
//                    }
//                }
//            }
//        }
//        d.neighborCount -= count;
        probe(d, newestSlideIndex, true);
    }

    private void processExpiredData(int expiredSlide) {

        all_slides.remove(expiredSlide);
        outlierList.remove(expiredSlide);

//        for (ArrayList<Vector> av : all_slides.values()) {
//            for (int i = 0; i < av.size(); i++) {
//                Vector d = av.get(i);
//                if (d.neighbor_count_map.containsKey(expiredSlide)) {
//                    d.neighborCount -= d.neighbor_count_map.get(expiredSlide);
//                    d.neighbor_count_map.remove(expiredSlide);
//                }
//                if (d.corePoints_map.containsKey(expiredSlide)) {
//                    d.corePoints_map.remove(expiredSlide);
//                }
//            }
//        }
        if (neighborCountTrigger.containsKey(expiredSlide)) {
            for (Vector d : neighborCountTrigger.get(expiredSlide)) {
                d.neighbor_count_map.remove(expiredSlide);
            }
        }
        if (neighborCountTrigger.containsKey(expiredSlide - 1)) {
            neighborCountTrigger.remove(expiredSlide - 1);
        }

        if (corePointMapTrigger.containsKey(expiredSlide)) {
            for (Vector d : corePointMapTrigger.get(expiredSlide)) {
                d.corePoints_map.remove(expiredSlide);
            }
        }
//        if(corePointMapTrigger.containsKey(expiredSlide))
        corePointMapTrigger.remove(expiredSlide - 1);

        all_core_points.remove(expiredSlide);

    }

    public double check_double_distance_neighbor(Vector d, Vector d2) {

        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R * 2) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R * 2) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R * 2) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R * 2) {
                return ub_d;
            }
        }

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_double_distance_neighbor_reversed(Vector d, Vector d2) {
//        if (Constants.useLB1) {
//            double lb_d = DistanceFunction.LB_distance1(d, d2);
//            if (lb_d > Constants.R) {
//                return lb_d;
//            }
//        } else if (Constants.useLB2) {
//            double lb_d = DistanceFunction.LB_distance2(d, d2);
//            if (lb_d > Constants.R) {
//                return lb_d;
//            }
//        }
//
//        if (Constants.useUB1) {
//            double ub_d = DistanceFunction.UB_distance1(d, d2);
//            if (ub_d <= Constants.R) {
//                return ub_d;
//            }
//        } else if (Constants.useUB2) {
//            double ub_d = DistanceFunction.UB_distance2(d, d2);
//            if (ub_d <= Constants.R) {
//                return ub_d;
//            }
//        }

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_R2neighbor(Vector d, Vector d2) {
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R / 2) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R / 2) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R / 2) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R / 2) {
                return ub_d;
            }
        }

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_neighbor(Vector d, Vector d2) {
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R) {
                return ub_d;
            }
        }

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_neighbor_reversed(Vector d, Vector d2) {

//        double ub_d = DistanceFunction.UB_distance2(d, d2);
//        if (ub_d <= Constants.R) {
//            return ub_d;
//        }
//
//        double lb_d = DistanceFunction.LB_distance2(d, d2);
//        if (lb_d > Constants.R) {
//            return lb_d;
//        }
        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

//    public double check_distance_neighbor2(Vector d, Vector d2) {
////        double lb_d = DistanceFunction.LB_distance2(d, d2);
////        if (lb_d > Constants.R) {
////            return lb_d;
////        }
////        double ub_d = DistanceFunction.UB_distance2(d, d2);
////        if (ub_d <= Constants.R) {
////            return ub_d;
////        }
//        double exact_d = DistanceFunction.euclideanDistance(d, d2);
//        return exact_d;
//    }
//    public void neighborUpdate(Vector d, int numNeighbor, int slideIndex) {
//        if (d.getSlideIndex() > slideIndex) {
//            if (d.preceding_neighbor_map.containsKey(slideIndex)) {
//                d.preceding_neighbor_map.put(slideIndex,
//                        d.preceding_neighbor_map.get(slideIndex) + numNeighbor);
//            } else {
//                d.preceding_neighbor_map.put(slideIndex, numNeighbor);
//                //add d to the trigger list
//                if (trigger_list.containsKey(slideIndex)) {
//                    trigger_list.get(slideIndex).add(d);
//                }
//            }
//        }
//    }
//
//    public void neighborUpdate(Vector d, List<Vector> neighbors) {
//        if (neighbors.isEmpty()) {
//            return;
//        }
//
//        Vector d2 = neighbors.get(0);
//        if (d.getSlideIndex() > d2.getSlideIndex()) {
//            if (d.preceding_neighbor_map.containsKey(d2.getSlideIndex())) {
//                d.preceding_neighbor_map.put(d2.getSlideIndex(),
//                        d.preceding_neighbor_map.get(d2.getSlideIndex()) + neighbors.size());
//            } else {
//                d.preceding_neighbor_map.put(d2.getSlideIndex(), neighbors.size());
//                //add d to the trigger list
//                if (trigger_list.containsKey(d2.getSlideIndex())) {
//                    trigger_list.get(d2.getSlideIndex()).add(d);
//                }
//            }
//
//        }
////        
//        for (Vector d3 : neighbors) {
//            if (d3.lastProb < d.getSlideIndex()
//                    && d3.getSlideIndex() < d.getSlideIndex()) {
//
//                if (d3.passiveNeighbor_map.containsKey(d.getSlideIndex())) {
//                    d3.passiveNeighbor_map.put(d.getSlideIndex(),
//                            d3.passiveNeighbor_map.get(d.getSlideIndex()) + 1);
//                } else {
//                    d3.passiveNeighbor_map.put(d.getSlideIndex(), 1);
//                }
////                d3.numNeighborPassive += 1;
//                d3.neighborCount += 1;
//            }
//        }
//    }
//    public void neighborUpdate(Vector d, Vector d2) {
//        if (d.getSlideIndex() > d2.getSlideIndex()) {
//            if (d.preceding_neighbor_map.containsKey(d2.getSlideIndex())) {
//                d.preceding_neighbor_map.put(d2.getSlideIndex(),
//                        d.preceding_neighbor_map.get(d2.getSlideIndex()) + 1);
//            } else {
//                d.preceding_neighbor_map.put(d2.getSlideIndex(), 1);
//            }
//        }
//
//        if (d2.lastProb > -1 && d2.lastProb < d.getSlideIndex()) {
//            d2.numNeighborPassive += 1;
//            d2.neighborCount += 1;
//
//        }
//    }
    public ArrayList<Vector> findNeighbors(Vector d, ArrayList<Vector> candidates,
            int maxReturn, int order) {
        ArrayList<Vector> results = new ArrayList<>();

        switch (order) {
            case 0:
                //using lb first
                for (Vector d2 : candidates) {
                    d2.checkedLB = false;
                    d2.checkedUB = false;
                }
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }
                    double lb_d = DistanceFunction.LB_distance2(d, d2);
                    if (lb_d > Constants.R) {
                        d2.checkedLB = true;

                    }
                }
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }
                    if (!d2.checkedLB
                            && DistanceFunction.UB_distance2(d, d2) <= Constants.R) {
                        results.add(d2);
                        d2.checkedUB = true;
                        if (results.size() >= maxReturn) {
                            return results;
                        }
                    }

                }
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }
                    if (!d2.checkedLB && !d2.checkedUB) {
                        if (DistanceFunction.euclideanDistance(d, d2) <= Constants.R) {
                            results.add(d2);
                            if (results.size() >= maxReturn) {
                                return results;
                            }
                        }
                    }
                }
                break;
            case 1:
                for (Vector d2 : candidates) {
                    d2.checkedLB = false;
                    d2.checkedUB = false;
                }
                //using ub first
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }
                    if (DistanceFunction.UB_distance2(d, d2) <= Constants.R) {
                        results.add(d2);
                        d2.checkedUB = true;
                        if (results.size() >= maxReturn) {
                            return results;
                        }
                    }

                }
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }

                    if (!d2.checkedUB && DistanceFunction.LB_distance2(d, d2) > Constants.R) {
                        d2.checkedLB = true;

                    }
                }
                for (Vector d2 : candidates) {
                    if (d2.arrivalTime == d.arrivalTime) {
                        continue;
                    }
                    if (d2.added_for_probing) {
                        continue;
                    }
                    if (!d2.checkedLB && !d2.checkedUB) {
                        if (DistanceFunction.euclideanDistance(d, d2) <= Constants.R) {
                            results.add(d2);
                            if (results.size() >= maxReturn) {
                                return results;
                            }
                        }
                    }
                }
                break;
            case 2:
                for (Vector v : candidates) {
                    if (!v.added_for_probing && v.arrivalTime != d.arrivalTime) {
                        double distance = check_distance_neighbor(d, v);
                        if (distance <= Constants.R) {
                            results.add(v);
                            if (results.size() >= maxReturn) {
                                return results;
                            }
                        }
                    }

                }
                break;
            default:
                break;
        }
        return results;

    }

    public void updateNeighborWithCore(Vector d, Vector c, int idx) {
        ArrayList<Vector> neighbors = c.closeNeighbors;
        if (idx == d.getSlideIndex()) {

            neighbors.remove(d);

        }
//        if (neighbors.size() > Constants.k) {
//            neighbors = new ArrayList<>(neighbors.subList(0, Constants.k));
//        }

//                        d.neighbor_map.put(idx, neighbors);
        d.neighbor_count_map.put(idx, neighbors.size());

        //put d to trigger list 
        if (neighborCountTrigger.containsKey(idx)) {
            neighborCountTrigger.get(idx).add(d);
        } else {
            HashSet<Vector> hs = new HashSet<>();
            hs.add(d);
            neighborCountTrigger.put(idx, hs);
        }
        if (idx >= d.getSlideIndex()) {
            d.numSucceedingNeighbors += neighbors.size();
        }
        d.neighborCount += neighbors.size();
        d.corePoints_map.put(idx, c);

        //add d to trigger list
        if (corePointMapTrigger.containsKey(idx)) {
            corePointMapTrigger.get(idx).add(d);
        } else {
            HashSet<Vector> hs = new HashSet<>();
            hs.add(d);
            corePointMapTrigger.put(idx, hs);
        }

        //update passive neighbor for close neighbor of c
//        if (idx < d.getSlideIndex()) {
//            for (Vector d2 : c.closeNeighbors) {
//                if (d2.lastProbCoreRight < d.getSlideIndex() && d2.lastProbRight < d.getSlideIndex()) {
//                    //add d to passive map
//                    if (d2.passiveNeighbor_map.containsKey(d.getSlideIndex())) {
//                        d2.passiveNeighbor_map.put(d.getSlideIndex(), d2.passiveNeighbor_map.get(d.getSlideIndex()) + 1);
//                        d2.neighborCount += 1;
//                    } else {
//                        d2.passiveNeighbor_map.put(d.getSlideIndex(), 1);
//                        d2.neighborCount += 1;
//                    }
//                }
//            }
//        }

    }

    public void probe_CorePoints(Vector d, int minSlideIndex, int maxSlideIndex, int direction, int maxDepth) {
        int idx;
        boolean condition;
        if (direction == 0) {//increasing 
            idx = minSlideIndex;
            condition = (d.isOutlier() && idx <= maxSlideIndex && idx <= minSlideIndex + maxDepth);
        } else {
            idx = maxSlideIndex;
            condition = (d.isOutlier() && idx >= minSlideIndex && idx >= maxSlideIndex - maxDepth);
        }

        while (condition) {
            if (d.passiveNeighbor_map.containsKey(idx)) {
                d.neighborCount -= d.passiveNeighbor_map.get(idx);
                d.passiveNeighbor_map.remove(idx);
            }
            if (d.corePoints_map.get(idx) != null) {
                Vector c = d.corePoints_map.get(idx);

                updateNeighborWithCore(d, c, idx);
            } else {

                ArrayList<Vector> core_points = all_core_points.get(idx);
                if (core_points != null) {
                    for (Vector c : core_points) {
//                    double distance = DistanceFunction.euclideanDistance(d, c);
                        double distance = check_distance_R2neighbor(d, c);
                        if (distance <= Constants.R / 2) {
                            updateNeighborWithCore(d, c, idx);
                            break;
                        }
                    }
                }
            }

            if (direction == 0) {
                d.lastProbCoreRight = idx;
            } else {
                d.lastProbCoreLeft = idx;
            }
            if (!d.isOutlier()) {
                break;
            }
            //update idx and condition
            if (direction == 0) {
                idx += 1;
                condition = (d.isOutlier() && idx <= maxSlideIndex && idx <= minSlideIndex + maxDepth);
            } else {
                idx -= 1;
                condition = (d.isOutlier() && idx >= minSlideIndex && idx >= maxSlideIndex - maxDepth);
            }

        }

    }

    public void probe(Vector d, int newestSlide, boolean reprobing) {
        countProb += 1;
        //probe core point first
        int direction = 1;
        if (!reprobing) {

            int maxSlide = d.lastProbCoreLeft - 1;
            if (d.lastProbCoreLeft == -1) {
                maxSlide = newestSlide;
            }

            int minSlide = Math.max(0, expiredSlideIndex + 1);
            direction = 1;
            probe_CorePoints(d, minSlide, maxSlide, direction, Constants.maxDepth);
        } else {
            //probe core points
            int minSlide = d.lastProbCoreRight + 1;
            direction = 0;
            int maxSlide = newestSlide;
            probe_CorePoints(d, minSlide, maxSlide, direction, Constants.maxDepth);
            if (d.isOutlier()) {
                direction = 1;
                maxSlide = d.lastProbCoreLeft - 1;
                minSlide = Math.max(0, expiredSlideIndex + 1);
                probe_CorePoints(d, minSlide, maxSlide, direction, Constants.maxDepth);
            }
        }

        if (!d.isOutlier()) {
            return;
        }

        //probe normal points
        int slideIndex = d.getSlideIndex();
        if (!reprobing) {
            if (d.lastProbLeft > -1) {
                slideIndex = d.lastProbLeft - 1;
            }
//            if (d.isCore && newestSlide == d.getSlideIndex()) {
//                slideIndex = d.getSlideIndex() - 1;
//            }

        } else if (reprobing) {
            if (d.lastProbRight != -1) {
                slideIndex = d.lastProbRight + 1;
            } else {
                slideIndex = d.getSlideIndex();
//                if (d.isCore) {
//                    slideIndex = d.getSlideIndex() + 1;
//                }
                if (d.lastProbLeft > -1 && d.lastProbLeft < d.getSlideIndex()) {
                    slideIndex = d.getSlideIndex() + 1;
                }
            }
        }

        boolean continueCondition = (d.isOutlier() && slideIndex > -1 && slideIndex > expiredSlideIndex);
        if (reprobing) {
            continueCondition = (d.isOutlier() && slideIndex <= newestSlide);
        }

        while (continueCondition) {

//            if (d.passiveNeighbor_map.containsKey(slideIndex)) {
//                d.neighborCount -= d.passiveNeighbor_map.get(slideIndex);
//                d.passiveNeighbor_map.remove(slideIndex);
//            }
            //check with core points first
            ArrayList<Vector> candidates = new ArrayList<>();
            boolean foundCore = false;
            Vector core = null;
//            double distance_to_core = 0;

            if (d.corePoints_map.containsKey(slideIndex)) {
                core = d.corePoints_map.get(slideIndex);
                foundCore = true;
                candidates = core.doubleRCandidates;
            } 
            else {
                ArrayList<Vector> core_points = all_core_points.get(slideIndex);

                if (core_points != null) {
                    for (int i = 0; i < core_points.size(); i++) {
                        Vector c = core_points.get(i);
//                    double distance = DistanceFunction.euclideanDistance(d, c);
                        double distance = check_distance_neighbor(d, c);
                        if (distance <= Constants.R) {
                            
                            
                            if (d.neighbor_count_map.containsKey(slideIndex)) {

                                d.neighbor_count_map.put(slideIndex,
                                        d.neighbor_count_map.get(slideIndex) + 1);

                            } else {

                                d.neighbor_count_map.put(slideIndex, 1);
                            }

                            //put d into trigger list 
                            if (neighborCountTrigger.containsKey(slideIndex)) {
                                neighborCountTrigger.get(slideIndex).add(d);
                            } else {
                                HashSet<Vector> hs = new HashSet<>();
                                hs.add(d);
                                neighborCountTrigger.put(slideIndex, hs);
                            }

                            if (slideIndex >= d.getSlideIndex()) {
                                d.numSucceedingNeighbors += 1;
                            }
                            d.neighborCount += 1;
                            candidates = c.doubleRCandidates;
                            core = c;
                            foundCore = true;
                            break;

                        }

                    }
                }
            }

            if (!foundCore) {
                candidates = all_slides.get(slideIndex);
            }
//            System.out.println("Candidate size = " + candidates.size());

            for (int i = 0; i < candidates.size(); i++) {
                candidates.get(i).added_for_probing = false;
            }
            //set all current neighbors in the same slide to addedtoprobing
            if (d.corePoints_map.containsKey(slideIndex)) {
                d.corePoints_map.get(slideIndex).added_for_probing = true;
                for (Vector v : d.corePoints_map.get(slideIndex).closeNeighbors) {
                    v.added_for_probing = true;
                }
            }
            
            //filter far core and its close neighbors
            for(Vector c: all_core_points.get(slideIndex)){
                if(DistanceFunction.euclideanDistance(d, c) > Constants.R*3/2){
                    c.added_for_probing = true;
                    for(Vector v: c.closeNeighbors){
                        v.added_for_probing = true;
                    }
                }
            }

            //binary search in distances if found core
//            if (foundCore) {
//                int index;
//                if (core.distance_to_neighbors.length >= Constants.k
//                        && core.distance_to_neighbors[Constants.k - 1] <= Constants.R - distance_to_core) {
//                    index = Constants.k;
//                } else if (core.distance_to_neighbors.length < Constants.k && core.distance_to_neighbors.length > 0
//                        && core.distance_to_neighbors[core.distance_to_neighbors.length - 1] <= Constants.R - distance_to_core) {
//                    index = core.distance_to_neighbors.length;
//                } else {
//
//                    index = Arrays.binarySearch(core.distance_to_neighbors, Constants.R - distance_to_core);
//                    if (index < 0) {
//                        index = -index - 1;
//
//                    }
//                }
////                System.out.println("Index = "+ index);
//                if (index >= Constants.k) {
//                    index = Constants.k;
//                } else if (index < 0) {
//                    index = 0;
//                }
//                d.neighborCount += index;
//                neighborUpdate(d, candidates.subList(0, index));
////                for (int i = 0; i < index; i++) {
////                    candidates.get(i).added_for_probing = true;
////
////                }
//                if (d.neighborCount >= Constants.k) {
//                    break;
//
//                } else {
//                    for (int i = index + 1; i < candidates.size(); i++) {
//                        Vector d2 = candidates.get(i);
//
//                        if (d2.arrivalTime != d.arrivalTime) {
////                            if (d2.added_for_probing) {
////                                continue;
////                            }
////                            d2.added_for_probing = true;
//                            double distance;
//                            distance = check_distance_neighbor(d, d2);
//
//                            if (distance <= Constants.R) {
//                                d.neighborCount += 1;
//                                neighborUpdate(d, d2);
//                                if (d.neighborCount >= Constants.k) {
//                                    break;
//
//                                }
//                            } else if (distance - distance_to_core >= Constants.R) {
//                                break;
//                            }
//                        }
//                    }
//                }
//
//            } else {
//            int countNeighborBefore = 0;
//            for (int s : d.neighbor_count_map.keySet()) {
//                if (s < slideIndex) {
//                    countNeighborBefore += d.neighbor_count_map.get(s);
//                }
//            }
            ArrayList<Vector> neighbors = findNeighbors(d, candidates,
                    Constants.k, 2);
            d.neighborCount += neighbors.size();

            if (slideIndex >= d.getSlideIndex()) {
                d.numSucceedingNeighbors += neighbors.size();
            }
            if (d.neighbor_count_map.containsKey(slideIndex)) {
                d.neighbor_count_map.put(slideIndex,
                        d.neighbor_count_map.get(slideIndex) + neighbors.size());
            } else {
                d.neighbor_count_map.put(slideIndex, neighbors.size());
            }

            //put d into trigger list 
            if (neighborCountTrigger.containsKey(slideIndex)) {
                neighborCountTrigger.get(slideIndex).add(d);
            } else {
                HashSet<Vector> hs = new HashSet<>();
                hs.add(d);
                neighborCountTrigger.put(slideIndex, hs);
            }
//            neighborUpdate(d, neighbors);
//            for (int i = 0; i < candidates.size(); i++) {
//                Vector d2 = candidates.get(i);
//
//                if (d2.arrivalTime != d.arrivalTime) {
//                    if (d2.added_for_probing) {
//                        continue;
//                    }
//                    d2.added_for_probing = true;
//                    double distance;
//                    distance = check_distance_neighbor(d, d2);
//
//                    if (distance <= Constants.R) {
//                        d.neighborCount += 1;
//                        neighborUpdate(d, d2);
//                        d2.added_for_probing = true;
//                        if (d.neighborCount >= Constants.k) {
//                            break;
//
//                        }
//                    }
////                        else if (distance - distance_to_core >= Constants.R) {
////                            break;
////                        }
//                }
//
//            }

            if (!d.isOutlier()) {
                break;

            }

            if (!reprobing) {

                d.lastProbLeft = slideIndex;
            } else {
                d.lastProbRight = slideIndex;
            }

            //update continuecondition
            if (!reprobing) {
                slideIndex--;
                continueCondition = (d.isOutlier() && slideIndex > -1 && slideIndex > expiredSlideIndex);
            } else {
                slideIndex++;
                continueCondition = (d.isOutlier() && slideIndex <= newestSlide);
            }
        }

        if (d.isOutlier()) {
            if (outlierList.containsKey(d.getSlideIndex())) {
                outlierList.get(d.getSlideIndex()).add(d);
            } else {
                HashSet hs = new HashSet<>();
                hs.add(d);
                outlierList.put(d.getSlideIndex(), hs);
            }
        } else if (d.numSucceedingNeighbors >= Constants.k) {
            d.neighbor_count_map.clear();
            d.corePoints_map.clear();
        }

        if (!reprobing) {
            avg_depth = (avg_depth * countProb + d.getSlideIndex() - slideIndex) / (countProb + 1);

        }

//        System.out.println("Count stored neighbors = "+ d.countStoredNeighbor());
    }

//    public static <T extends Number> int[] asArray(final T... a) {
//        int[] b = new int[a.length];
//        for (int i = 0; i < b.length; i++) {
//            b[i] = a[i].intValue();
//        }
//        return b;
//    }
//
//    public static int[] argsort(final ArrayList<Double> a, final boolean ascending) {
//        Integer[] indexes = new Integer[a.size()];
//        for (int i = 0; i < indexes.length; i++) {
//            indexes[i] = i;
//        }
//        Arrays.sort(indexes, new Comparator<Integer>() {
//            @Override
//            public int compare(final Integer i1, final Integer i2) {
//                return (ascending ? 1 : -1) * Double.compare(a.get(i1), a.get(i2));
//            }
//        });
//        return asArray(indexes);
//    }
    public ArrayList<Vector> selectCore(int slideIndex) {
        ArrayList<Vector> corePoints = new ArrayList<>();
        ArrayList<Vector> possibleCandidates = all_slides.get(slideIndex);

        for (Vector c : possibleCandidates) {
            if (!c.checked) {
                c.checked = true;
                c.isCore = true;
                c.lastProbLeft = c.getSlideIndex();
                c.lastProbCoreRight = c.getSlideIndex();
                c.lastProbRight = c.getSlideIndex();
                c.lastProbCoreLeft = c.getSlideIndex();

                corePoints.add(c);
                for (Vector d : possibleCandidates) {
                    if (c.arrivalTime != d.arrivalTime) {
                        double distance = check_double_distance_neighbor(d, c);
                        if (distance <= Constants.R * 2) {
                            c.doubleRCandidates.add(d);
//                    c.distance_to_neighbors.add(distance);

                            if (distance <= Constants.R) {
                                c.neighborCount += 1;
                                d.checked = true;

                                if (!d.passiveNeighbor_map.containsKey(slideIndex)) {
                                    d.passiveNeighbor_map.put(slideIndex, 1);
                                } else {
                                    d.passiveNeighbor_map.put(slideIndex, 1 + d.passiveNeighbor_map.get(slideIndex));
                                }
                                d.neighborCount += 1;

                                if (distance <= Constants.R / 2) {
//                            if(c.closeNeighbors.size() < Constants.k)
                                    c.closeNeighbors.add(d);
                                    //add c to corePointMap of d
                                    d.corePoints_map.put(slideIndex, c);

//                                    if (c.closeNeighbors.size() >= Constants.k) {
//                                        c.isSafe = true;
//                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
        return corePoints;
    }

}
