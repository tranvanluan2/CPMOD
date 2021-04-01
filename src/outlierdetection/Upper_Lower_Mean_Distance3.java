/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class Upper_Lower_Mean_Distance3 {

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public static HashMap<Integer, ArrayList<CorePoint>> all_core_points = new HashMap<>();

    public static HashMap<Integer, HashSet<C_Data>> outlierList = new HashMap<>();

    public static HashMap<Integer, HashSet<CorePoint>> corePointMapTrigger = new HashMap<>();
    public static HashMap<Integer, HashSet<C_Data>> neighborCountTrigger = new HashMap<>();

    public double timeFindingCore = 0;
    public double timeCheckingCandidates = 0;
    public double timeAddingToNeighborCount = 0;
    public double timeMarkProbing = 0;
    public double avg_candidate_size = 0;

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        currentTime = _currentTime;
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W / slide;
        System.out.println("Expire slide index = " + expiredSlideIndex);

        processExpiredData(expiredSlideIndex);

        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
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

        for (int i = 0; i < data.size(); i++) {
            Data o = data.get(i);
            C_Data d = new C_Data(o);

            d_to_process.add(d);

            ArrayList<C_Data> idx = all_slides.get(d.sIndex);
            if (idx != null) {

                idx.add(d);
            } else {
                idx = new ArrayList<>(Constants.slide);
                idx.add(d);
                all_slides.put(d.sIndex, idx);
            }
        }

        for (Integer sIdx : slide_to_process) {

            ArrayList<CorePoint> corePoints = selectCore(sIdx);
            all_core_points.put(sIdx, corePoints);
            System.out.println("Num new core = " + corePoints.size());

        }

        System.out.println("Probing");
        long start = Utils.getCPUTime();
        int newestSlide = (currentTime - 1) / Constants.slide;
        timeCheckingCandidates = 0;
        timeFindingCore = 0;
        timeAddingToNeighborCount = 0;
        timeMarkProbing = 0;
        for (int i = 0; i < d_to_process.size(); i++) {
            C_Data d = d_to_process.get(i);
            if (d.neighborCount < Constants.k) {
                probe(d, newestSlide);
            }
//            System.out.println("Finished "+ i);
        }
        System.out.println("Time finding cores = " + timeFindingCore);
        System.out.println("Time checking candidates = " + timeCheckingCandidates);
        System.out.println("Time adding to neighbor counts = " + timeAddingToNeighborCount);
        System.out.println("Time mark probing = " + timeMarkProbing);
        System.out.println("Time for probing = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        System.out.println("re-probing + outlier report");
//        start = Utils.getCPUTime();
        for (int slideIndex : outlierList.keySet()) {
            for (C_Data d : outlierList.get(slideIndex)) {
                if (d.neighborCount < Constants.k && d.sIndex < newestSlide) {
                    if (d.lastProbRight < newestSlide) {
                        probe(d, newestSlide);
                    }

                }
                if (d.neighborCount < Constants.k) {
                    result.add(d);
                }
            }
        }
//        System.out.println("Time for reprobing+report outliers = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        return result;

    }

    public boolean check_distance_neighbor_boolean(C_Data d, C_Data d2) {
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R) {
                return false;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R) {
                return false;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R) {
                return true;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R) {
                return true;
            }
        }

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d <= Constants.R;
    }

    public double check_distance_neighbor(C_Data d, C_Data d2) {
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

    private void processExpiredData(int expiredSlideIndex) {
        all_slides.remove(expiredSlideIndex);
        outlierList.remove(expiredSlideIndex);
        if (neighborCountTrigger.containsKey(expiredSlideIndex)) {
            for (C_Data d : neighborCountTrigger.get(expiredSlideIndex)) {
                d.neighborCount -= d.pred_neighbor_count.get(expiredSlideIndex);
                d.pred_neighbor_count.remove(expiredSlideIndex);
            }
        }
        if (neighborCountTrigger.containsKey(expiredSlideIndex - 1)) {
            neighborCountTrigger.remove(expiredSlideIndex - 1);
        }

        all_core_points.remove(expiredSlideIndex);
    }

    private ArrayList<CorePoint> selectCore(Integer sIdx) {

        ArrayList<CorePoint> corePoints = new ArrayList<>();
        ArrayList<C_Data> possibleCandidates = all_slides.get(sIdx);

        for (int i = 0; i < Constants.slide; i++) {
            C_Data d = all_slides.get(sIdx).get(i);
            Random r = new Random();
            if ((d.closeCoreMaps_R.isEmpty() && d.closeCoreMaps_halfR.isEmpty())) {

                CorePoint c = new CorePoint(d);

                all_slides.get(sIdx).set(i, c);
                corePoints.add(c);
                c.lastProbLeft = c.sIndex;
                c.lastProbRight = c.sIndex;
                for (C_Data d2 : possibleCandidates) {
                    if (d2.arrivalTime != c.arrivalTime) {
                        double distance = DistanceFunction.euclideanDistance(c, d2);

                        if (distance <= Constants.R / 2) {
                            c.neighborCount += 1;
                            c.closeNeighbors_halfR.add(d2);
                            if (d2.closeCoreMaps_halfR.containsKey(sIdx)) {
                                d2.closeCoreMaps_halfR.get(sIdx).add(c);
                            } else {
                                ArrayList<CorePoint> corepoints = new ArrayList<>();
                                corepoints.add(c);
                                d2.closeCoreMaps_halfR.put(sIdx, corepoints);
                            }
                        } else if (distance <= Constants.R) {
                            c.closeNeighbors_R.add(d2);
                            c.neighborCount += 1;
//                                    d.neighborCount += 1;
                            if (d2.closeCoreMaps_R.containsKey(sIdx)) {
                                d2.closeCoreMaps_R.get(sIdx).add(c);
                            } else {
                                ArrayList<CorePoint> corepoints = new ArrayList<>();
                                corepoints.add(c);
                                d2.closeCoreMaps_R.put(sIdx, corepoints);
                            }
                        } else if (distance <= Constants.R * 1.5) {
                            c.closeNeighbors_3halfR.add(d2);
                            if (d2.closeCoreMaps_3halfR.containsKey(sIdx)) {
                                d2.closeCoreMaps_3halfR.get(sIdx).add(c);
                            } else {
                                ArrayList<CorePoint> corepoints = new ArrayList<>();
                                corepoints.add(c);
                                d2.closeCoreMaps_3halfR.put(sIdx, corepoints);
                            }
                        } else if (distance <= Constants.R * 2) {
                            c.closeNeighbors_2R.add(d2);
                            if (d2.closeCoreMaps_2R.containsKey(sIdx)) {
                                d2.closeCoreMaps_2R.get(sIdx).add(c);
                            } else {
                                ArrayList<CorePoint> corepoints = new ArrayList<>();
                                corepoints.add(c);
                                d2.closeCoreMaps_2R.put(sIdx, corepoints);
                            }
                        }

                    }

                }

            }
        }
        return corePoints;

    }

    private void probe_slide(C_Data d, int slideIndex) {

        int countNeighbor = 0;

        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();

        //find close core
        long start = Utils.getCPUTime();
        ResultFindCore rf = findCloseCore(d, slideIndex);
//        System.out.println("Time to find core = "+ (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        timeFindingCore += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        if (cores != null) {
            if (distance <= Constants.R / 2) {
                CorePoint c = cores.get(0);

                //grab close neighbor in range R/2 of c
                countNeighbor += c.closeNeighbors_halfR.size();
                if (countNeighbor >= Constants.k 
                        || d.numSucceedingNeighbor + c.closeNeighbors_halfR.size() >= Constants.k) {
                    d.neighborCount += c.closeNeighbors_halfR.size();
                    if (slideIndex < d.sIndex) {
                        

                        d.pred_neighbor_count.put(slideIndex, c.closeNeighbors_halfR.size());
                        if (neighborCountTrigger.containsKey(slideIndex)) {
                            neighborCountTrigger.get(slideIndex).add(d);
                        } else {
                            HashSet<C_Data> hs = new HashSet<>();
                            hs.add(d);
                            neighborCountTrigger.put(slideIndex, hs);
                        }
                    } else {
                        d.numSucceedingNeighbor += c.closeNeighbors_halfR.size();
                    }

                    return;

                }
                possibleCandidates.add(c.closeNeighbors_R);
                possibleCandidates.add(c.closeNeighbors_3halfR);

            } else if (distance <= Constants.R) {
                Collections.sort(cores, new Comparator<CorePoint>() {
                    @Override
                    public int compare(CorePoint o1, CorePoint o2) {
                        if (o1.getTotalCoverPoint() < o2.getTotalCoverPoint()) {
                            return -1;
                        } else if (o1.getTotalCoverPoint() == o2.getTotalCoverPoint()) {
                            return 0;
                        } else {
                            return 1;
                        }
                    }

                });

                possibleCandidates.add(cores.get(0).closeNeighbors_halfR);
                possibleCandidates.add(cores.get(0).closeNeighbors_R);
                possibleCandidates.add(cores.get(0).closeNeighbors_3halfR);
                possibleCandidates.add(cores.get(0).closeNeighbors_2R);

            } else if (distance <= Constants.R * 2) {
                for (CorePoint c : cores) {
                    possibleCandidates.add(c.closeNeighbors_R);
                    possibleCandidates.add(c.closeNeighbors_halfR);

                    possibleCandidates.add(c.closeNeighbors_3halfR);
                    possibleCandidates.add(c.closeNeighbors_2R);
                }
            }

        } 
//        else {
////            System.out.println("AAAAAAAAAAAAAAAAa");
////            possibleCandidates.add(all_slides.get(slideIndex));
//        }

        //check with other points 
        //first mark the probed points
//        start = Utils.getCPUTime();
//        for (C_Data d2 : possibleCandidates) {
//            d2.addedForProbing = false;
//        }
//
//        if (cores != null && distance <= Constants.R / 2) {
//            for (C_Data d2 : cores.get(0).closeNeighbors_halfR) {
//                d2.addedForProbing = true;
//            }
//        }
//        timeMarkProbing += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//        System.out.println("Time for marking addedForProbing = "+ timeMarkProbing);
        //filtering first
//        ArrayList<CorePoint> corePoints = all_core_points.get(slideIndex);
//        for (CorePoint c2 : corePoints) {
//            distance = DistanceFunction.euclideanDistance(d, c2);
//            if (distance > Constants.R * 3) {
//                for (C_Data d2 : c2.closeNeighbors_2R) {
//                    d2.addedForProbing = true;
//                }
//            } 
//            else if (distance > Constants.R * 2.5) {
//                for (C_Data d2 : c2.closeNeighbors_3halfR) {
//                    d2.addedForProbing = true;
//                }
//            } else if (distance > Constants.R * 2) {
//                for (C_Data d2 : c2.closeNeighbors_R) {
//                    d2.addedForProbing = true;
//                }
//            } else if (distance > Constants.R * 1.5) {
//                for (C_Data d2 : c2.closeNeighbors_halfR) {
//                    d2.addedForProbing = true;
//                }
//
//            }
//        }
        start = Utils.getCPUTime();
        for (ArrayList<C_Data> ps : possibleCandidates) {
            for (C_Data d2 : ps) {
//                if (!d2.addedForProbing) {
//                distance = DistanceFunction.euclideanDistance(d, d2);

                if (check_distance_neighbor_boolean(d, d2)) {
                    countNeighbor += 1;
                    if(d2.sIndex >= d.sIndex){
                        d.numSucceedingNeighbor +=1;
                    }
                    if (countNeighbor >= Constants.k || d.numSucceedingNeighbor >= Constants.k) {

                        break;
                    }
                }
//                    d2.addedForProbing = true;
//                }
            }
            if (countNeighbor >= Constants.k) {

                break;
            }
        }
        timeCheckingCandidates += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        d.neighborCount += countNeighbor;
        if (slideIndex < d.sIndex) {
//            d.numSucceedingNeighbor += countNeighbor;

            d.pred_neighbor_count.put(slideIndex, countNeighbor);
//            start = Utils.getCPUTime();
            if (neighborCountTrigger.containsKey(slideIndex)) {
                neighborCountTrigger.get(slideIndex).add(d);
            } else {
                HashSet<C_Data> hs = new HashSet<>();
                hs.add(d);
                neighborCountTrigger.put(slideIndex, hs);
            }
//            timeAddingToNeighborCount += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        }
//        System.out.println("Time looping through candidates = "+ ((Utils.getCPUTime() - start) * 1.0 / 1000000000));

    }

    private void probe(C_Data d, int newestSlide) {
        if (d.lastProbRight < newestSlide) {
            //prob right first
            int slideIndex = d.lastProbRight + 1;
            if (d.lastProbRight == -1) {
                slideIndex = d.sIndex;
            }
            while (slideIndex <= newestSlide && d.neighborCount < Constants.k) {
                probe_slide(d, slideIndex);
                d.lastProbRight = slideIndex;
                slideIndex++;
            }
        }
        //prob left
        if (d.neighborCount < Constants.k) {
            int slideIndex = d.lastProbLeft - 1;
            if (d.lastProbLeft == -1) {
                slideIndex = d.sIndex - 1;
            }
            while (slideIndex > expiredSlideIndex && slideIndex >= 0
                    && d.neighborCount < Constants.k) {
                probe_slide(d, slideIndex);
                d.lastProbLeft = slideIndex;
                slideIndex--;
            }
        }

        if (d.neighborCount < Constants.k) {
            //add to outlier List
            if (outlierList.containsKey(d.sIndex)) {
                outlierList.get(d.sIndex).add(d);
            } else {
                HashSet hs = new HashSet();
                hs.add(d);
                outlierList.put(d.sIndex, hs);
            }
        }
    }

    final class ResultFindCore {

        private final double distance;
        private final ArrayList<CorePoint> cores;

        public ResultFindCore(double distance, ArrayList<CorePoint> cores) {
            this.distance = distance;
            this.cores = cores;
        }

        public double getDistance() {
            return this.distance;
        }

        public ArrayList<CorePoint> getCore() {
            return this.cores;
        }
    }

    private ResultFindCore findCloseCore(C_Data d, int slideIndex) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR.containsKey(slideIndex)) {
            resultCore = d.closeCoreMaps_halfR.get(slideIndex);
            return new ResultFindCore(Constants.R / 2, resultCore);
        } else if (d.closeCoreMaps_R.containsKey(slideIndex)) {
            resultCore = d.closeCoreMaps_R.get(slideIndex);
            return new ResultFindCore(Constants.R, resultCore);
        }

        ArrayList<CorePoint> corePoints = all_core_points.get(slideIndex);

        ArrayList<CorePoint> inRangeRCores = new ArrayList<>();
        ArrayList<CorePoint> inRangeDoubleRCores = new ArrayList<>();
        for (CorePoint c : corePoints) {
            double distance = DistanceFunction.euclideanDistance(d, c);
            if (distance <= Constants.R / 2) {
                resultCore = new ArrayList<>();
                resultCore.add(c);
                return new ResultFindCore(Constants.R / 2, resultCore);
            } else if (distance <= Constants.R) {
                inRangeRCores.add(c);
                //test
//                return new ResultFindCore(Constants.R, inRangeRCores);
                //end test
            } else if (distance <= Constants.R * 2) {
                inRangeDoubleRCores.add(c);
            }
        }
        if (!inRangeRCores.isEmpty()) {
//            System.out.println("in Range R core = "+ inRangeRCores.size());

            return new ResultFindCore(Constants.R, inRangeRCores);
        } else if (!inRangeDoubleRCores.isEmpty()) {
            return new ResultFindCore(Constants.R * 2, inRangeDoubleRCores);
        } else {
            return new ResultFindCore(Constants.R * 2, null);
        }
    }

    class C_Data extends Data {

//        public boolean addedForProbing = false;
        private int numSucceedingNeighbor = 0;
//        private boolean isOutlier;
        public int lastProbRight = -1;
        public int lastProbLeft = -1;
//        private int lastProbCoreLeft = -1;
//        private int lastProbCoreRight = -1;
        private HashMap<Integer, Integer> pred_neighbor_count = new HashMap<>();
        public int neighborCount = 0;

        private HashMap<Integer, ArrayList<CorePoint>> closeCoreMaps_halfR = new HashMap<>();
        private HashMap<Integer, ArrayList<CorePoint>> closeCoreMaps_R = new HashMap<>();
        private HashMap<Integer, ArrayList<CorePoint>> closeCoreMaps_3halfR = new HashMap<>();
        private HashMap<Integer, ArrayList<CorePoint>> closeCoreMaps_2R = new HashMap<>();
//        private CorePoint closeCore;
//        public ArrayList<CorePoint> linkCoreMaps = new ArrayList<>();

        public ArrayList<Integer[]> groups = new ArrayList<>();

        public double[] mean_list;
        public double[] std_list;
        public ArrayList<Integer> ordered_attributes = new ArrayList<>();

        public int sIndex = -1;

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;

            if (null != Constants.dataFile) {
                switch (Constants.dataFile) {

                    case "household2.txt":

                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6});
                        }
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(2, 3, 6, 0, 5, 1, 4));
                        }
                        break;
                    case "new_hpc.txt":
                        this.groups.add(new Integer[]{0, 1, 3});
                        this.groups.add(new Integer[]{2, 4, 5, 6});
                        break;
                    case "covtype.data":
                        if (Constants.useLB1 || Constants.useUB1) {
//                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
//                            this.groups.add(new Integer[]{ });
                            this.groups.add(new Integer[]{0});
                            this.groups.add(new Integer[]{1});
                            this.groups.add(new Integer[]{2});
                            this.groups.add(new Integer[]{4});
                            this.groups.add(new Integer[]{3});
                            this.groups.add(new Integer[]{5});
                            this.groups.add(new Integer[]{6});
                            this.groups.add(new Integer[]{7});
                            this.groups.add(new Integer[]{8});
                            this.groups.add(new Integer[]{9});
                            this.groups.add(new Integer[]{11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
//                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5});
//                            this.groups.add(new Integer[]{6, 7, 8, 9, 54});
                        }

                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54));
                        }
                        break;
                    case "ethylene.txt":

                        if (Constants.useLB1 || Constants.useUB1) {

                            this.groups.add(new Integer[]{0, 2});
//                    this.groups.add(new Integer[]{1});
                            this.groups.add(new Integer[]{3, 4, 6, 7});
                            this.groups.add(new Integer[]{8, 14, 15});
                            this.groups.add(new Integer[]{1, 9});
                            this.groups.add(new Integer[]{10, 11, 12});
                            this.groups.add(new Integer[]{13});
                            this.groups.add(new Integer[]{5});
                        }

                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7));
                        }
                        break;

                    case "new_tao.txt":
                        this.groups.add(new Integer[]{0, 1, 2});

                        break;

                    case "tao.txt":
                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[]{0, 1, 2});
                        }
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(1, 2));
                        }
                        break;

                    case "gaussian.txt":
                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[0]);
                        }
                    default:
                        break;
                }

                if (this.groups.size() > 0) {
                    int num_group = this.groups.size();

                    this.mean_list = new double[num_group];
                    this.std_list = new double[num_group];
                    for (int i = 0; i < num_group; i++) {
                        double temp = 0;
                        for (Integer j : this.groups.get(i)) {
                            temp += this.values[j];
                        }
                        this.mean_list[i] = temp / this.groups.get(i).length;
                        if (this.groups.get(i).length == 1) {
                            this.std_list[i] = 0;
                        } else {
                            temp = 0;
                            for (Integer j : this.groups.get(i)) {
                                temp += (this.values[j] - this.mean_list[i]) * (this.values[j] - this.mean_list[i]);
                            }

                            this.std_list[i] = Math.sqrt(temp / this.groups.get(i).length);
                        }
                    }
                }
            }

//            this.isOutlier = true;
        }

        public C_Data() {

        }

        public int countNeighbor() {
            return neighborCount;
        }

    }

    class CorePoint extends C_Data {

        public ArrayList<C_Data> closeNeighbors_halfR = new ArrayList<>();
        public ArrayList<C_Data> closeNeighbors_R = new ArrayList<>();
        public ArrayList<C_Data> closeNeighbors_3halfR = new ArrayList<>();
        public ArrayList<C_Data> closeNeighbors_2R = new ArrayList<>();

        public int getTotalCoverPoint() {
            return closeNeighbors_2R.size() + closeNeighbors_3halfR.size()
                    + closeNeighbors_R.size() + closeNeighbors_halfR.size();
        }
//  
//        private int numPointsCovered = 0;

        public CorePoint(C_Data d) {
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.arrivalTime = d.arrivalTime;
            this.groups = d.groups;
            this.mean_list = d.mean_list;
            this.std_list = d.std_list;
        }

    }

}
