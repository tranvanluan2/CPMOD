/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import be.tarsos.lsh.Vector;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author Luan Tran
 */
public class Upper_Lower_Mean_Distance_2 {

//    public static PriorityQueue<Vector> event_queue = new PriorityQueue(new VectorNeighborComparator());
    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<Vector>> all_slides = new HashMap<>();
//    public static ArrayList<Integer> number_points_in_clusters = new ArrayList<>();

//    public static HashMap<Integer, ArrayList<Vector>> all_core_points = new HashMap<>();
    public static PriorityQueue<Vector> all_core_points = new PriorityQueue<>(new Comparator<Vector>() {
        @Override
        public int compare(Vector o1, Vector o2) {
            if (o1.arrivalTime > o2.arrivalTime) {
                return -1;
            } else {
                return 1;
            }
        }
    });

//    public static HashMap<Integer, HashSet<Vector>> trigger_list = new HashMap<>();
//    public static HashMap<Integer, HashSet<Vector>> outlierList = new HashMap<>();
//    public static int countCheckNeighbor = 0;
//    public static int countFiltedByLB = 0;
////    public static int countFiltedByUB = 0;
//
//    public static double timeForProcessingExpiredSlide = 0;
//    public static double timeForProcessingNewData = 0;
//    public static double timeForProcessingReprobing = 0;
//    public static int numFoundCore = 0;
    public Upper_Lower_Mean_Distance_2() {

    }

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
//        numFoundCore = 0;
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        currentTime = _currentTime;
        expiredSlideIndex = (currentTime - Constants.W - 1) / slide;

        long start = Utils.getCPUTime();
        processExpiredData(expiredSlideIndex);
        System.out.println("Time for processing expired data = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

//        start = Utils.getCPUTime();
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
            Random r = new Random();

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
        if (data.size() == Constants.W) {
            selectCore();
        } else {
            start = Utils.getCPUTime();
            updateCore(d_to_process.get(0).getSlideIndex());
            System.out.println("Time for updating core = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        }

        System.out.println("Number of cores = " + all_core_points.size());

        processExpiredCore(expiredSlideIndex);
        int newestSlide = (currentTime - 1) / slide;

        for (Integer sIdx : all_slides.keySet()) {
            for (Vector v : all_slides.get(sIdx)) {
                if (v.isCore) {
                    v.numTotalCloseNeighbor = v.totalCloseNeighbor();
                }
            }
        }
        start = Utils.getCPUTime();
        for (Integer sIdx : all_slides.keySet()) {
            for (Vector d : all_slides.get(sIdx)) {
                if (d.isCore) {
                    if (d.neighborCount < Constants.k) {
                        result.add(d);
                    }
                } else {
                    if (d.closeToFullCore()) {
                        continue;
                    }
                    if (d.numSucceedingNeighbors < Constants.k
                            && d.totalNeighborCount() < Constants.k) {
                        probeForNeighbors(d, newestSlide);
                    }
                    if (d.numSucceedingNeighbors < Constants.k
                            && d.totalNeighborCount() < Constants.k) {
                        result.add(d);
                    }
                }
            }
        }
        System.out.println("Time for probing = " + ((Utils.getCPUTime() - start) * 1.0 / 1000000000));

        return result;

    }

    private void processExpiredCore(int expiredSlide) {
        System.out.println("Processing expired core....");
        if (expiredSlide > -1) {
            ArrayList<Vector> to_remove = new ArrayList<>();

//            for (int i = all_core_points.size() - 1; i >= 0; i--) {
            for (Vector c : all_core_points) {

//                c.numCloseNeighborMap.remove(expiredSlide);
                c.largeRNeighborCandidates.remove(expiredSlide);
                c.linkedPointMaps.remove(expiredSlide);
                if (c.getSlideIndex() <= expiredSlide && c.linkedPointMaps.isEmpty()) {
                    to_remove.add(c);
                }

            }
//            }
            for (Vector c : to_remove) {

                c.largeRNeighborCandidates.clear();
//                c.numCloseNeighborMap.clear();

                all_core_points.remove(c);
            }
        }
    }

    private void processExpiredData(int expiredSlide) {
        System.out.println("Processing expired data....");
        ArrayList<Vector> expireData = all_slides.get(expiredSlide);
        if (expireData != null) {
            for (Vector d : expireData) {
                d.closeNeighborCores.clear();
                d.neighbor_count_map.clear();
//                d.closeNeighborCores = null;

            }
        }
        all_slides.remove(expiredSlide);
        //remove expired data from core point neighbors/candidates

        if (expiredSlide > -1) {

            for (Integer sIdx : all_slides.keySet()) {
                for (Vector v : all_slides.get(sIdx)) {
                    v.neighbor_count_map.remove(expiredSlide);

                }
            }
        }

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

    public void probe_slide(Vector d, int sIdx) {

        int countNeighborSoFar = 0;
        //select core with the smallest size 
        Vector bestCore = null;
        int countCandidates = Integer.MAX_VALUE;
        for (Vector v : d.closeNeighborCores) {
            int count = 0;

            if (v.numCloseNeighborMap.containsKey(sIdx)) {
                count -= v.numCloseNeighborMap.get(sIdx);
            }

            if (count < countCandidates) {
                countCandidates = count;
                bestCore = v;
            }
        }

        //probe using the candidates from v
//        if (d.closeNeighborCores == bestCore) {
        if (bestCore.linkedPointMaps.containsKey(sIdx)) {
            if (d.getSlideIndex() == sIdx) {
                d.neighbor_count_map.put(sIdx, bestCore.linkedPointMaps.get(sIdx).size() - 1);
                countNeighborSoFar = bestCore.linkedPointMaps.get(sIdx).size() - 1;
            } else {
                d.neighbor_count_map.put(sIdx, bestCore.linkedPointMaps.get(sIdx).size());
                countNeighborSoFar = bestCore.linkedPointMaps.get(sIdx).size();
            }
        }

//        }
        if (countNeighborSoFar >= Constants.k) {
            return;
        }

        if (sIdx >= d.getSlideIndex()) {
            d.numSucceedingNeighbors += countNeighborSoFar;
            if (d.numSucceedingNeighbors >= Constants.k) {

                return;
            }
        }
        if (bestCore.largeRNeighborCandidates.containsKey(sIdx)) {
            for (Vector d2 : bestCore.largeRNeighborCandidates.get(sIdx)) {
                if (d2.closeNeighborCores.contains(bestCore)) {
                    continue;
                }
//                double distance = DistanceFunction.euclideanDistance(d, d2);
                double distance = check_distance_neighbor(d, d2);
                if (distance <= Constants.R) {
                    if (d.neighbor_count_map.containsKey(sIdx)) {
                        d.neighbor_count_map.put(sIdx, d.neighbor_count_map.get(sIdx) + 1);
                    } else {
                        d.neighbor_count_map.put(sIdx, 1);
                    }
                    countNeighborSoFar++;
                    if (countNeighborSoFar >= Constants.k) {
                        return;
                    }

                    if (sIdx >= d.getSlideIndex()) {
                        d.numSucceedingNeighbors += 1;
                        if (d.numSucceedingNeighbors >= Constants.k) {
                            //d is safe inlier
//                            d.neighborCores.clear();
//                            d.closeNeighborCores = null;
//                            d.neighbor_count_map.clear();
                            return;
                        }
                    }
                }
            }
        }

    }

    public void probeForNeighbors(Vector d, int newestSlide) {
        if (d.lastProbRight <= newestSlide) {
            //can Prob Right

            int sIdx = d.lastProbRight + 1;
            if (d.lastProbRight == -1) {
                sIdx = d.getSlideIndex();
            }
            while (sIdx <= newestSlide) {

                probe_slide(d, sIdx);

                d.lastProbRight = sIdx;

                sIdx++;

                if (d.totalNeighborCount() > Constants.k) {
                    break;
                }
            }
        }

        if (d.totalNeighborCount() < Constants.k && d.lastProbLeft == -1) {
            int sIdx = d.getSlideIndex() - 1;
            while (sIdx > expiredSlideIndex && sIdx >= 0) {
                probe_slide(d, sIdx);
                d.lastProbLeft = sIdx;
                sIdx--;
                if (d.totalNeighborCount() > Constants.k) {
                    break;
                }

            }
        }
    }

    public boolean coreNeedUpdate(Vector c) {

        if (c.getSlideIndex() > expiredSlideIndex) {

            if (c.numSucceedingNeighbors < Constants.k) {
                return true;
            }
        }

//        boolean hasUnSafeInlier = false;
//        for(Integer sidx: c.linkedPointMaps.keySet()){
//            for(Vector d: c.linkedPointMaps.get(sidx)){
//                if(d.numSucceedingNeighbors < Constants.k) {
//                    hasUnSafeInlier = true;
//                    break;
//                }
//                
//            }
//            if(hasUnSafeInlier) break;
//        }
        for (Integer sidx : c.linkedPointMaps.keySet()) {

            if (c.linkedPointMaps.get(sidx).size() < Constants.k) {
                return true;
            }
        }
        return false;
    }

    public void updateCore(int slideIndex) {
        //to be 
        ArrayList<Vector> new_data = all_slides.get(slideIndex);
        long start = Utils.getCPUTime();
        List<Integer> all_slide_idxs = new ArrayList<>(all_slides.keySet());
        Collections.sort(all_slide_idxs, Collections.reverseOrder());

        int countNewCore = 0;

        for (Vector d : new_data) {
            for (Vector c : all_core_points) {
                if (!c.linkedPointMaps.containsKey(slideIndex)
                        || c.linkedPointMaps.get(slideIndex).size() < Constants.k
                        || !d.neighbor_count_map.containsKey(slideIndex)
                        || d.neighbor_count_map.get(slideIndex) < Constants.k) {
                    update_point_with_core(d, c);
                }

            }
        }
        for (Vector v : new_data) {
            if (!v.isCore && v.closeNeighborCores.isEmpty()) {
                //promote d to be a core
                v.isCore = true;
                boolean quit = false;
                for (Integer sIndex : all_slide_idxs) {
                    int countCloseNeighbor = 0;
                    for (Vector d : all_slides.get(sIndex)) {
                        if (v.arrivalTime != d.arrivalTime
                                && (!v.linkedPointMaps.containsKey(slideIndex)
                                || v.linkedPointMaps.get(slideIndex).size() < Constants.k
                                || d.numSucceedingNeighbors < Constants.k)) {
                            countCloseNeighbor += update_point_with_core(d, v);

                        } else if (v.linkedPointMaps.containsKey(slideIndex)
                                && v.linkedPointMaps.get(slideIndex).size() >= Constants.k
                                && d.numSucceedingNeighbors >= Constants.k) {
                            quit = true;
                            break;
                        }
                    }
                    if (quit) {
                        break;
                    }

                }
            }
        }

    }

    public int update_point_with_core(Vector d, Vector v) {
        int result = 0;
        int sIndex = d.getSlideIndex();
//        double distance = check_double_distance_neighbor(d, v);
        double distance = DistanceFunction.euclideanDistance(d, v);
        if (distance <= Constants.R * 3 / 2) {
            if (v.largeRNeighborCandidates.containsKey(sIndex)) {
                v.largeRNeighborCandidates.get(sIndex).add(d);
            } else {
                ArrayList<Vector> candidateList = new ArrayList<>();
                candidateList.add(d);
                v.largeRNeighborCandidates.put(sIndex, candidateList);
            }

            if (distance <= Constants.R) {
                v.neighborCount += 1;
                if (d.getSlideIndex() >= v.getSlideIndex()) {
                    v.numSucceedingNeighbors += 1;
                }
                if (v.neighbor_count_map.containsKey(d.getSlideIndex())) {
                    v.neighbor_count_map.put(d.getSlideIndex(), v.neighbor_count_map.get(d.getSlideIndex()) + 1);
                } else {
                    v.neighbor_count_map.put(d.getSlideIndex(), 1);
                }

                if (distance <= Constants.R / 2) {
                    d.closeNeighborCores.add(v);

                    if (v.linkedPointMaps.containsKey(d.getSlideIndex())) {
                        v.linkedPointMaps.get(d.getSlideIndex()).add(d);
                    } else {
                        ArrayList<Vector> links = new ArrayList<>();
                        links.add(d);
                        v.linkedPointMaps.put(d.getSlideIndex(), links);
                    }

                    result = 1;

                }
            }

        }
        return result;
    }

    public void selectCore() {

        List<Integer> all_slide_idxs = new ArrayList<>(all_slides.keySet());
        Collections.sort(all_slide_idxs, Collections.reverseOrder());
        for (Integer slideIndex : all_slide_idxs) {
            for (Vector v : all_slides.get(slideIndex)) {
                //
                if (v.closeNeighborCores.isEmpty()) {
                    v.isCore = true;
                    all_core_points.add(v);
                    for (Integer sIndex : all_slide_idxs) {

                        for (Vector d : all_slides.get(sIndex)) {
                            if (v.arrivalTime != d.arrivalTime) {
                                update_point_with_core(d, v);
//                                if (countCloseNeighbor >= Constants.k) {
//                                    break;
//                                }
                            }
                        }
//                        if (v.neighborCount >= Constants.k) {
//                            break;
//                        }

                    }
                }

            }

        }

    }

}
