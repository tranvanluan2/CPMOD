/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class New_MQ_CPOD_2 {

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();

//    public static ArrayList<CorePoint> all_distinct_core_map = new ArrayList<>();
    public static HashMap<Integer, HashMap<Double, ArrayList<CorePoint>>> all_core_points_map = new HashMap<>();
    public static ArrayList<CorePoint> all_distince_core_point = new ArrayList<>();
    public static ArrayList<Double> all_r_queries = new ArrayList<>();
//    public static double min_r = Double.MAX_VALUE;

    public ArrayList<Double> all_r;
    public ArrayList<Integer> all_k;
    public double timeForCreatCore = 0;
    public double timeForFindNeighbors = 0;
    public double timeForProcessingExpiredSlide = 0;
    public double count = 0;

    public HashMap<OD_Query, HashSet<Data>> slide_process(ArrayList<Data> data, int _currentTime) {

        HashMap<OD_Query, HashSet<Data>> result = new HashMap<>();
        for (OD_Query q : all_queries) {
            result.put(q, new HashSet<>());
        }
        currentTime = _currentTime;
        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
        int[] slide_to_process;
        if (currentTime == Constants.W) {
            slide_to_process = new int[Constants.W / Constants.slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
        } else {
            slide_to_process = new int[]{(currentTime - 1) / Constants.slide};
        }

        //insert new slide
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

        count += 1;

        double start = Utils.getCPUTime();

        expiredSlideIndex = (currentTime - 1) / Constants.slide - Constants.W / Constants.slide;

        processExpiredData(expiredSlideIndex);
        timeForProcessingExpiredSlide += (Utils.getCPUTime() - start) * 1.0 / 1000000000;

        start = Utils.getCPUTime();
        for (int sIdx : slide_to_process) {
            selectAllCore(sIdx, all_r_queries);
        }
        System.out.println("#distinct core points = " + all_distince_core_point.size());

        int newestSlide = (currentTime - 1) / Constants.slide;
        if (currentTime == Constants.W) {
            for (CorePoint c : all_distince_core_point) {
                for (Double r : c.supported_r) {
                    c.computeTotalHalfRPoints(r);
                }
            }
        } else {
            for (Double r : all_r_queries) {

                for (CorePoint c : all_core_points_map.get(newestSlide).get(r)) {

                    if (c.totalHalfRPoints.get(r) != null) {
                        c.updateHalfRPoints(newestSlide, r);
                    } else {
                        c.computeTotalHalfRPoints(r);
                    }

                }
            }
        }
        timeForCreatCore += (Utils.getCPUTime() - start) * 1.0 / 1000000000;

//        System.out.println("Print for testingm, for R = 1.9 ===============");
//
//        for (CorePoint c : all_core_points_map.get(0).get(1.9)) {
//            System.out.println(c.arrivalTime);
//            System.out.println("#half r points = " + c.totalHalfRPoints.get(1.9));
//        }
//        System.out.println("==========================");
//        System.out.println("Core Point for R = 0.19");
//        for (CorePoint c : all_core_points_map.get(0).get(0.19)) {
//            System.out.println(c.arrivalTime);
//        }
//        System.out.println("End testing ====================");
        start = Utils.getCPUTime();
        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                ArrayList<OD_Query> queries = d.getNeedProbQueries();

                if (!queries.isEmpty()) {
                    probe(d, newestSlide, queries);
                }
                for (OD_Query q : queries) {
                    if (result.containsKey(q)) {
                        result.get(q).add(d);
                    } else {
                        HashSet<Data> outliers = new HashSet<>();
                        outliers.add(d);
                        result.put(q, outliers);
                    }
                }
            }
        }
        timeForFindNeighbors += (Utils.getCPUTime() - start) * 1.0 / 1000000000;

        return result;

    }

    public void get_unique_k() {
        all_k = new ArrayList<>();
        for (OD_Query q : all_queries) {
            boolean exist = false;
            for (Integer k : all_k) {
                if (q.k == k) {
                    exist = true;
                    break;
                }
            }
            if (!exist) {
                all_k.add(q.k);
            }
        }
        Collections.sort(all_k);
    }

    public void get_unique_r() {
        all_r = new ArrayList<>();
        for (OD_Query q : all_queries) {
            boolean exist = false;
            for (Double r : all_r) {
                if (q.R == r) {
                    exist = true;
                    break;
                }
            }
            if (!exist) {
                all_r.add(q.R);
            }
        }
        Collections.sort(all_r);
    }

    public void add_query(OD_Query q) {
        all_queries.add(q);

        if (!all_r_queries.contains(q.R)) {
            all_r_queries.add(q.R);
        }

        //sort queries in increasing R order
        Collections.sort(all_queries, (OD_Query o1, OD_Query o2) -> {

            if (o1.R > o2.R) {
                return 2;
            } else if (o1.R == o2.R) {
                if (o1.k < o2.k) {
                    return 1;
                } else if (o1.k > o2.k) {
                    return -1;
                } else {
                    return 0;
                }
            } else {
                return -2;
            }
        });

        get_unique_r();
        get_unique_k();
    }

    public void remove_query(OD_Query q) {
        all_queries.remove(q);
    }

    public void remove_query(int idx) {
        all_queries.remove(idx);
    }

    public void update_query(int idx, OD_Query new_q) {
        all_queries.set(idx, new_q);
    }

    private double get_threshold(int W, int S, int m) {
        double delta = 12 * S * (W / m - 1);
        double k = -27 * 4 * W * S / (2 * Math.sqrt(Math.pow(delta, 3)));
        if (Math.abs(k) <= 1) {
            double x1 = 2 * Math.sqrt(delta) * Math.cos(Math.acos(k) / 3.0) / 3;
            double x2 = 2 * Math.sqrt(delta) * Math.cos(Math.acos(k) / 3.0 - 2 * Math.PI / 3) / 3;
            double x3 = 2 * Math.sqrt(delta) * Math.cos(Math.acos(k) / 3.0 + 2 * Math.PI / 3) / 3;

            return Math.max(x3, Math.max(x1, x2));
        } else {
            double x = Math.sqrt(delta) * Math.abs(k) / (3 * k) * (Math.pow(Math.abs(k) + Math.sqrt(k * k - 1), 1.0 / 3)
                    + Math.pow(Math.abs(k) - Math.sqrt(k * k - 1), 1.0 / 3));
            return x;
        }
    }

    private ArrayList<Double> get_threshold_new(int W, int S, double alpha, double beta, double gamma) {
        Cubic cubic = new Cubic();
        double a = 1;
        double b = 1;
        double c = S * alpha - W * gamma;
        double d = W * S * beta;
        cubic.solve(a, b, c, d);
        if (cubic.nRoots == 3) {
            return new ArrayList<>(Arrays.asList(cubic.x1, cubic.x2, cubic.x3));
//            return Math.max(Math.max(cubic.x1, cubic.x2), cubic.x3);
        } else {
            return new ArrayList<>(Arrays.asList(cubic.x1));
        }
    }

    private double get_dynamic_threshold(double r2, double r1, int W, int S, double T2r2, int m, int n) {
        double h = W * T2r2;
        double t = r1 * S / (r2 - r1);
        double a = 1;
        double b = 3 * t / m;
        double c = S - 3 * t + W * t / m - h;
        double d = W * S - W * t;
        Cubic cubic = new Cubic();
        cubic.solve(a, b, c, d);
        if (cubic.nRoots == 3) {
            return Math.max(Math.max(cubic.x1, cubic.x2), cubic.x3);
        } else {
            return cubic.x1;
        }
    }

    private void probe(C_Data d, int newestSlide, ArrayList<OD_Query> queries) {

        if (d.lastProbRight < newestSlide) {
            //prob right first
            int slideIndex = d.lastProbRight + 1;
            if (d.lastProbRight == -1) {
                slideIndex = d.sIndex;
            }

            while (slideIndex <= newestSlide && !queries.isEmpty()) {
                probe_slide_right2(d, slideIndex, queries);
                d.lastProbRight = slideIndex;
                slideIndex++;
            }

        }

        int slideIndex = d.lastProbLeft - 1;
        if (d.lastProbLeft == -1) {
            slideIndex = d.sIndex - 1;
        }

        while (slideIndex > expiredSlideIndex && slideIndex >= 0 && !queries.isEmpty()) {

            probe_slide_left2(d, slideIndex, queries);
            //probe left
            d.lastProbLeft = slideIndex;
            slideIndex--;
        }

        d.refineNeighborStore();
    }

    private void probe_slide_left2(C_Data d, int slideIndex, ArrayList<OD_Query> queries) {
        if (!d.pred_neighbor_count.containsKey(slideIndex)) {
            d.pred_neighbor_count.put(slideIndex, new HashMap<>());
        }

        HashMap<Double, Integer> pred_neigh = d.pred_neighbor_count.get(slideIndex);
        ArrayList<Double> all_supported_r = new ArrayList<Double>(all_core_points_map.get(slideIndex).keySet());

        Collections.sort(all_supported_r);

        ArrayList<Double> all_r_in_queries = new ArrayList<>();
        for (OD_Query q : queries) {
            if (!all_r_in_queries.contains(q.R)) {
                all_r_in_queries.add(q.R);
            }
        }

        Collections.sort(all_r_in_queries);

        Double max_r_in_query = all_r_in_queries.get(all_r_in_queries.size() - 1);

        boolean[] checked = new boolean[Constants.slide];

        int count_checked = 0;
        int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

        for (Double target_r : all_supported_r) {
//        while (!queries.isEmpty() && !all_r_in_queries.isEmpty()
//                && count_checked < Constants.slide) {
            if (queries.isEmpty() || count_checked >= Constants.slide || all_r_in_queries.isEmpty()) {
                break;
            }
            //scan possible points 
            ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();

            //find close core
            ResultFindCore rf = findCloseCore(d, slideIndex, target_r);
            double distance = rf.getDistance();
            ArrayList<CorePoint> cores = rf.getCore();
            if (cores != null) {
                if (distance <= target_r / 2) {
                    CorePoint c = cores.get(0);
                    possibleCandidates.add(c.getDataInRange(0, 1.5 * target_r, slideIndex));

                } else if (distance <= target_r) {

                    possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * target_r, slideIndex));

                } else if (distance <= target_r * 2) {
                    for (CorePoint c : cores) {
                        possibleCandidates.add(c.getDataInRange(0, target_r, slideIndex));

                    }

                }

            }

            for (ArrayList<C_Data> ps : possibleCandidates) {

                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
                    if (!checked[d2.arrivalTime - min_arrival_time]) {
                        double dist = DistanceFunction.euclideanDistance(d, d2);
                        for (Double r : all_r_in_queries) {

                            if (dist <= r) {
                                if (pred_neigh.containsKey(r) && pred_neigh.get(r) + get_neighbor_count(d.numSucceedingNeighbor, r)
                                        < all_k.get(all_k.size() - 1)) {
                                    pred_neigh.put(r, pred_neigh.get(r) + 1);
                                } else {
                                    pred_neigh.put(r, 1);
                                }

                                //check neighbor count, remove queries if d is inlier
                                ArrayList<OD_Query> to_remove = new ArrayList<>();
                                for (OD_Query q : queries) {
                                    if (q.R == r) {
                                        if (pred_neigh.get(r)
                                                + get_neighbor_count(d.numSucceedingNeighbor, r) >= q.k) {
                                            to_remove.add(q);
                                        }
                                    }
                                }
                                for (OD_Query q : to_remove) {
                                    queries.remove(q);
                                }

                            }
                        }

                        checked[d2.arrivalTime - min_arrival_time] = true;
                        count_checked += 1;
                    }
                    
                    if(queries.isEmpty() || all_r_in_queries.isEmpty()) break;
                }

            }
            //remove  r
            while(!all_r_in_queries.isEmpty() && all_r_in_queries.get(0) <= target_r){
                all_r_in_queries.remove(0);
            }
           
            if (target_r >= max_r_in_query) {
                break;
            }
        }
    }

    private void probe_slide_left(C_Data d, int slideIndex, ArrayList<OD_Query> queries, double max_r) {
//        HashSet<Double> all_r_in_queries = new HashSet<>();
//        for (OD_Query q : queries) {
//            all_r_in_queries.add(q.R);
//        }

        ArrayList<Double> all_r_in_queries = new ArrayList<>();
        for (OD_Query q : queries) {
            if (!all_r_in_queries.contains(q.R)) {
                all_r_in_queries.add(q.R);
            }
        }

        Collections.sort(all_r_in_queries);
        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();
        //find close core
        ResultFindCore rf = findCloseCore(d, slideIndex, max_r);
        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= max_r / 2) {
                CorePoint c = cores.get(0);
                possibleCandidates.add(c.getDataInRange(0, 1.5 * max_r, slideIndex));

            } else if (distance <= max_r) {

                possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * max_r, slideIndex));

            } else if (distance <= max_r * 2) {
                case_ = 1;

                for (CorePoint c : cores) {
                    possibleCandidates.add(c.getDataInRange(0, max_r, slideIndex));

                }

            }

            int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

            boolean[] checked = null;
            if (case_ == 1) {
                checked = new boolean[Constants.slide];
            }

            HashMap<Double, Integer> oldNeighborCount = d.neighborCounts();
//            for (Double r : all_r_in_queries) {
//                if (d.neighborCount.containsKey(r)) {
//                    oldNeighborCount.put(r, d.neighborCount.get(r));
//                } else {
//                    oldNeighborCount.put(r, 0);
//                }
//            }
//        outerloop:
            HashMap<Double, Integer> newNeighborCount = d.neighborCounts();
            for (ArrayList<C_Data> ps : possibleCandidates) {
                if (newNeighborCount.get(all_r.get(0)) >= all_k.get(all_k.size() - 1)) {
                    break;
                }
                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
//                if (!checked.contains(d2)) {
                    if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {
                        double dist = DistanceFunction.euclideanDistance(d, d2);
                        for (Double r : all_r_in_queries) {
                            if (dist <= r && newNeighborCount.get(r) < all_k.get(all_k.size() - 1)
                                    + get_neighbor_count(d.numSucceedingNeighbor, r)) {

                                if (newNeighborCount.containsKey(r)) {
                                    newNeighborCount.put(r, newNeighborCount.get(r) + 1);
                                } else {
                                    newNeighborCount.put(r, 1);
                                }

                            }
                        }

                        if (case_ == 1) {
                            checked[d2.arrivalTime - min_arrival_time] = true;
                        }
                    }
                }

            }
            for (Double r : all_r_in_queries) {
                if (!d.pred_neighbor_count.containsKey(slideIndex)) {
                    d.pred_neighbor_count.put(slideIndex, new HashMap<>());
                }
                HashMap<Double, Integer> pred_neigh = d.pred_neighbor_count.get(slideIndex);
                if (pred_neigh.containsKey(r)) {
                    pred_neigh.put(r, pred_neigh.get(r) + newNeighborCount.get(r) - oldNeighborCount.get(r));
                } else {
                    pred_neigh.put(r, newNeighborCount.get(r) - oldNeighborCount.get(r));
                }
                //check neighbor count, remove queries if d is inlier
                ArrayList<OD_Query> to_remove = new ArrayList<>();
                for (OD_Query q : queries) {
                    if (q.R == r) {
                        if (newNeighborCount.get(r) >= q.k) {
                            to_remove.add(q);
                        }
                    }
                }
                for (OD_Query q : to_remove) {
                    queries.remove(q);
                }
            }
        } else {
//            System.out.println("No core found!!!");
        }
    }

    private void processExpiredData(int expiredSlideIndex) {
//        if (all_slides.containsKey(expiredSlideIndex)) {
//            for (C_Data d : all_slides.get(expiredSlideIndex)) {
//                d.reset();
//            }
//        }
        all_slides.remove(expiredSlideIndex);
        for (CorePoint c : all_distince_core_point) {
            for (Double r : c.supported_r) {
                c.discountHalfPoints(expiredSlideIndex, r);
            }
            for (Bin b : c.all_bins) {
                if (b.data.containsKey(expiredSlideIndex)) {
                    b.data.remove(expiredSlideIndex);
                }
            }

        }

        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                if (d.pred_neighbor_count.containsKey(expiredSlideIndex)) {
//                    HashMap<Double, Integer> pred_count = d.pred_neighbor_count.get(expiredSlideIndex);
                    //update neighbor count
//                    for (Double r : all_r) {
//                        Integer pred = get_neighbor_count(pred_count, r);
//                        d.neighborCount.put(r, d.neighborCount.get(r) - pred);
//                    }

                    d.pred_neighbor_count.remove(expiredSlideIndex);
                }
            }
        }
        if (all_core_points_map.containsKey(expiredSlideIndex)) {
            all_core_points_map.get(expiredSlideIndex).clear();
        }
        all_core_points_map.remove(expiredSlideIndex);

        ArrayList<CorePoint> to_remove = new ArrayList<>();
        for (CorePoint c : all_distince_core_point) {
            if (c.supported_slides <= expiredSlideIndex) {
                to_remove.add(c);
            }
        }
        for (CorePoint c : to_remove) {
            c.reset();
            all_distince_core_point.remove(c);
        }
    }

    private Integer get_neighbor_count(HashMap<Double, Integer> pred_count, Double r) {
        if (pred_count.containsKey(r)) {
            return pred_count.get(r);
        }
        //find closest smaller r
        for (int i = all_r.size() - 1; i >= 0; i--) {
            if (all_r.get(i) < r && pred_count.containsKey(all_r.get(i))) {
                return pred_count.get(all_r.get(i));
            }
        }

        return 0;
    }

    final class ResultFindCore {

        private final double distance;
        private final ArrayList<CorePoint> cores;
        public ArrayList<Double> distance_to_cores = new ArrayList<>();

        public ResultFindCore(double distance, ArrayList<CorePoint> cores) {
            this.distance = distance;
            this.cores = cores;
        }

        public ResultFindCore(double distance, ArrayList<CorePoint> cores, ArrayList<Double> all_distances) {
            this.distance = distance;
            this.cores = cores;
            this.distance_to_cores = all_distances;
        }

        public double getDistance() {
            return this.distance;
        }

        public ArrayList<CorePoint> getCore() {
            return this.cores;
        }
    }

    private ResultFindCore findCloseCore(C_Data d, int slideIndex, double radius) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR.get(radius) != null
                && all_core_points_map.get(slideIndex).get(radius).contains(d.closeCoreMaps_halfR.get(radius))) {
            resultCore = new ArrayList<>();
            resultCore.add(d.closeCoreMaps_halfR.get(radius));
            return new ResultFindCore(DistanceFunction.euclideanDistance(d, d.closeCoreMaps_halfR.get(radius)), resultCore);
        } else if (d.closeCoreMaps_R.get(radius) != null) {
            CorePoint c = d.closeCoreMaps_R.get(radius);
            if (all_core_points_map.get(slideIndex).get(radius).contains(c)) {
                resultCore = new ArrayList<>();
                resultCore.add(c);
                return new ResultFindCore(DistanceFunction.euclideanDistance(d, c), resultCore);
            }

        }

        ArrayList<CorePoint> corePoints = all_core_points_map.get(slideIndex).get(radius);

        ArrayList<CorePoint> inRangeRCores = new ArrayList<>();
        ArrayList<CorePoint> inRangeDoubleRCores = new ArrayList<>();
        ArrayList<Double> distance_to_cores = new ArrayList<>();

        if (corePoints != null) {
            for (int i = 0; i < corePoints.size(); i++) {
                CorePoint c = corePoints.get(i);
                double distance = DistanceFunction.euclideanDistance(d, c);

                if (distance <= radius) {
                    inRangeRCores.add(c);

                    break;
                } else if (distance <= radius * 2) {
                    inRangeDoubleRCores.add(c);
                    distance_to_cores.add(distance);
                }
            }
        }
        if (!inRangeRCores.isEmpty()) {
            return new ResultFindCore(radius, inRangeRCores);
        } else if (!inRangeDoubleRCores.isEmpty()) {

            return new ResultFindCore(radius * 2, inRangeDoubleRCores, distance_to_cores);
        } else {
            return new ResultFindCore(radius * 2, null);
        }
    }

    private void probe_slide_right2(C_Data d, int slideIndex, ArrayList<OD_Query> queries) {

        ArrayList<Double> all_supported_r = new ArrayList<Double>(all_core_points_map.get(slideIndex).keySet());

        Collections.sort(all_supported_r);

        ArrayList<Double> all_r_in_queries = new ArrayList<>();
        for (OD_Query q : queries) {
            if (!all_r_in_queries.contains(q.R)) {
                all_r_in_queries.add(q.R);
            }
        }

        Collections.sort(all_r_in_queries);

        for (Double r : all_r_in_queries) {
            d.numSucceedingNeighbor.put(r, get_neighbor_count(d.numSucceedingNeighbor, r));
        }

        boolean[] checked = new boolean[Constants.slide];

        int count_checked = 0;
        int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

        Double max_r_in_query = all_r_in_queries.get(all_r_in_queries.size() - 1);

        for (Double target_r : all_supported_r) {
            if (queries.isEmpty() || count_checked >= Constants.slide || all_r_in_queries.isEmpty()) {
                break;
            }
            //scan possible points 
            ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();

            //find close core
            ResultFindCore rf = findCloseCore(d, slideIndex, target_r);
            double distance = rf.getDistance();
            ArrayList<CorePoint> cores = rf.getCore();
            if (cores != null) {
                if (distance <= target_r / 2) {
                    CorePoint c = cores.get(0);
                    possibleCandidates.add(c.getDataInRange(0, 1.5 * target_r, slideIndex));

                } else if (distance <= target_r) {

                    possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * target_r, slideIndex));

                } else if (distance <= target_r * 2) {
                    for (CorePoint c : cores) {
                        possibleCandidates.add(c.getDataInRange(0, target_r, slideIndex));

                    }

                }

            }

            for (ArrayList<C_Data> ps : possibleCandidates) {
                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
                    if (!checked[d2.arrivalTime - min_arrival_time]) {
                        double dist = DistanceFunction.euclideanDistance(d, d2);
                        for (Double r : all_r_in_queries) {
                            if (dist <= r) {
                                d.numSucceedingNeighbor.put(r, d.numSucceedingNeighbor.get(r) + 1);
                                ArrayList<OD_Query> to_remove = new ArrayList<>();
                                for (OD_Query q : queries) {
                                    if (q.R == r) {
                                        if (d.numSucceedingNeighbor.get(r) >= q.k) {
                                            d.safeInlierQueries.add(q);
                                            to_remove.add(q);
                                        }

                                    }
                                }
                                for (OD_Query q : to_remove) {
                                    queries.remove(q);
                                }

                            }
                        }

                        checked[d2.arrivalTime - min_arrival_time] = true;
                        count_checked += 1;
                        
                        
                    }
                    if(queries.isEmpty() || all_r_in_queries.isEmpty()) break;
                }

            }
            HashMap<Double, Integer> neighborCount = d.neighborCounts();
            ArrayList<OD_Query> to_remove = new ArrayList<>();
            for (OD_Query q : queries) {
                if (neighborCount.get(q.R) >= q.k) {
                    to_remove.add(q);
                }
            }
            for (OD_Query q : to_remove) {
                queries.remove(q);
            }
            //remove  r
            while(!all_r_in_queries.isEmpty() && all_r_in_queries.get(0) <= target_r){
                all_r_in_queries.remove(0);
            }
            
            if (target_r >= max_r_in_query) {
                break;
            }
        }

    }

    private void probe_slide_right(C_Data d, int slideIndex, ArrayList<OD_Query> queries, double max_r) {

        ArrayList<Double> all_r_in_queries = new ArrayList<>();
        for (OD_Query q : queries) {
            if (!all_r_in_queries.contains(q.R)) {
                all_r_in_queries.add(q.R);
            }
        }

        Collections.sort(all_r_in_queries);

//        HashSet<Double> all_r = new HashSet<>();
//        for (OD_Query q : queries) {
//            all_r.add(q.R);
//        }
        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();
        //find close core
//        long start = Utils.getCPUTime();
        ResultFindCore rf = findCloseCore(d, slideIndex, max_r);
        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= max_r / 2) {
                CorePoint c = cores.get(0);
                //grab close neighbor in range R/2 of c

                possibleCandidates.add(c.getDataInRange(0, 1.5 * max_r, slideIndex));

            } else if (distance <= max_r) {

                possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * max_r, slideIndex));

            } else if (distance <= max_r * 2) {
                case_ = 1;

                for (CorePoint c : cores) {
                    possibleCandidates.add(c.getDataInRange(0, max_r, slideIndex));

                }

            }

            int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

            boolean[] checked = null;
            if (case_ == 1) {
                checked = new boolean[Constants.slide];
            }

            HashMap<Double, Integer> oldNumSucNeighbor = new HashMap<>();
            for (Double r : all_r_in_queries) {
                if (d.numSucceedingNeighbor.containsKey(r)) {
                    oldNumSucNeighbor.put(r, d.numSucceedingNeighbor.get(r));
                } else {
                    oldNumSucNeighbor.put(r, 0);
                }
            }
//        outerloop:
            for (ArrayList<C_Data> ps : possibleCandidates) {

                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
//                if (!checked.contains(d2)) {
                    if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {
                        double dist = DistanceFunction.euclideanDistance(d, d2);
                        for (Double r : all_r) {
                            if (dist <= r) {
                                if (d.numSucceedingNeighbor.containsKey(r)) {
                                    d.numSucceedingNeighbor.put(r, d.numSucceedingNeighbor.get(r) + 1);
                                } else {
                                    d.numSucceedingNeighbor.put(r, 1);
                                }

                                ArrayList<OD_Query> to_remove = new ArrayList<>();
                                for (OD_Query q : queries) {
                                    if (q.R == r) {
                                        if (d.numSucceedingNeighbor.get(r) >= q.k) {
                                            d.safeInlierQueries.add(q);
                                            to_remove.add(q);
                                        }
                                    }
                                }
                                for (OD_Query q : to_remove) {
                                    queries.remove(q);
                                }

                            }
                        }

                        if (case_ == 1) {
                            checked[d2.arrivalTime - min_arrival_time] = true;
                        }
                    }
                }

            }

            HashMap<Double, Integer> neighborCount = d.neighborCounts();
            for (Double r : all_r_in_queries) {
//                neighborCount.put(r, neighborCount.get(r) + get_neighbor_count(d.numSucceedingNeighbor, r) - oldNumSucNeighbor.get(r));
                //check neighbor count, remove queries if d is inlier
                ArrayList<OD_Query> to_remove = new ArrayList<>();
                for (OD_Query q : queries) {
                    if (q.R == r) {
                        if (neighborCount.get(r) >= q.k) {
                            to_remove.add(q);
                        }
                    }
                }
                for (OD_Query q : to_remove) {
                    queries.remove(q);
                }
            }
        } else {
//            System.out.println("No core found!!!");
        }
    }

    private void selectAllCore(int sIdx, ArrayList<Double> all_r) {
        ArrayList<C_Data> data = all_slides.get(sIdx);
        //sort r in decreasing order
        Collections.sort(all_r);
        all_core_points_map.put(sIdx, new HashMap<>());
        for (Double r : all_r) {
            HashMap<Double, ArrayList<CorePoint>> map = all_core_points_map.get(sIdx);
            map.put(r, new ArrayList<>());
        }

        ArrayList<CorePoint> newCores = new ArrayList<>();

//        Double min_r = all_r.get(0);
//        Double max_r = all_r.get(all_r.size() - 1);
        for (C_Data d : data) {
            //find largest R that d is not connected to 
            //check with current core point 
            //find largest r the data point is not supported

            ArrayList<Double> r_not_found = new ArrayList<>();
            ArrayList<Double> r_found = new ArrayList<>();
            ArrayList<CorePoint> core_found = new ArrayList<>();
            ArrayList<Double> dists = new ArrayList<>();
            ArrayList<Boolean> isAdded = new ArrayList<>();
            for (int i = 0; i < all_r.size(); i++) {

                Double r = all_r.get(i);
                boolean found = false;
                for (CorePoint c : all_core_points_map.get(sIdx).get(r)) {

                    double dist = DistanceFunction.euclideanDistance(d, c);
                    if (dist <= r) {
                        found = true;
                        r_found.add(r);
                        core_found.add(c);
                        dists.add(dist);
                        isAdded.add(true);
                        break;
                    }

                }

                //find core in the history 
                if (!found) {
                    for (CorePoint c : all_distince_core_point) {
                        if (c.supported_r.contains(r)) {
                            double dist = DistanceFunction.euclideanDistance(d, c);

                            if (dist <= r) {
                                found = false;
                                r_found.add(r);
                                core_found.add(c);
                                dists.add(dist);
                                isAdded.add(false);
                                break;
                            }
                        }
                    }
                }
                if (!found) {
                    r_not_found.add(r);
                }

            }

            //create corepoint if r_not_found is not empty
            if (!r_not_found.isEmpty()) {
                //create core for r for that no core is found
                CorePoint c = new CorePoint(d);
                newCores.add(c);
                //connect d 
                c.creatBins(all_r);
//                all_distince_core_point.add(c);
                for (int i = 0; i < r_not_found.size(); i++) {
                    d.closeCoreMaps_halfR.put(r_not_found.get(i), c);

                    if (sIdx > c.supported_slides) {
                        c.supported_slides = sIdx;
                    }
                    all_core_points_map.get(sIdx).get(r_not_found.get(i)).add(c);

                    c.supported_r.add(r_not_found.get(i));
                }
                c.putDataToBin(d, 0, sIdx);

            }
            //connect to core point if r_found is not empty
            if (!r_found.isEmpty()) {
                for (int i = 0; i < r_found.size(); i++) {
                    Double r = r_found.get(i);
                    Double distance = dists.get(i);
                    boolean added = isAdded.get(i);
                    CorePoint core = core_found.get(i);
                    if (distance <= r / 2) {
                        d.closeCoreMaps_halfR.put(r, core);
                    } else {
                        d.closeCoreMaps_R.put(r, core);
                    }
                    if (!added) {
                        if (sIdx > core.supported_slides) {
                            core.supported_slides = sIdx;
                        }
                        all_core_points_map.get(sIdx).get(r).add(core);
                    }
                    if (!core_found.subList(0, i).contains(core)) {
                        core.putDataToBin(d, distance, sIdx);
                    }
                }
            }

        }

        //refine, remove core points for too small r
        HashSet<Double> r_to_remove = new HashSet<>();
//        for (int i = all_r.size() - 2; i >= 0; i--) {
//            Double r = all_r.get(i);
//            Double closest_r = -1.0;
//            for (int j = i + 1; j < all_r.size(); j++) {
//                if (all_core_points_map.get(sIdx).get(all_r.get(j)).size() > 0) {
//                    closest_r = all_r.get(j);
//                    break;
//                }
//            }
//            double threshold = get_threshold(Constants.W, Constants.slide, all_core_points_map.get(sIdx).get(closest_r).size());
////            System.out.println("Computing for r = " + r);
////            System.out.println("Threshold = "+ threshold);
//            
//            if (all_core_points_map.get(sIdx).get(r).size() > threshold) {
//                r_to_remove.add(r);
//                all_core_points_map.get(sIdx).get(r).clear();
////                System.out.println("Not creating cores for r = "+ r);
//            }
//            else{
////                System.out.println("++Keep creating cores for r = "+ r);
//            }
//        }
//        for (CorePoint c : newCores) {
//
//            for (Double r : r_to_remove) {
//                c.supported_r.remove(r);
//            }
//            if (!c.supported_r.isEmpty()) {
//                all_distince_core_point.add(c);
//            }
//
//        }
//        ArrayList<CorePoint> to_remove = new ArrayList<>();
//        for(CorePoint c: all_distince_core_point){
//            if(c.supported_r.isEmpty()){
//                to_remove.add(c);
//            }
//        }
//        for(CorePoint c: to_remove){
//            all_distince_core_point.remove(c);
//        }
        //
        double lastgamma = -1;
        Double lastR = null;
        double lastalpha = -1;
        double lastbeta = -1;
        boolean[] checked = new boolean[Constants.slide];
        HashSet<CorePoint> checked_core = new HashSet<>();
        for (int r_idx = all_r.size() - 1; r_idx >= 0; r_idx--) {
            Double radius = all_r.get(r_idx);
//        for (Double radius : all_r) {
            boolean needCore = true;
            if (lastgamma > 0) {
                //compute threshold 
                System.out.println("Last alpha = " + lastalpha);
                System.out.println("Last beta = " + lastbeta);
                System.out.println("Last gamma = " + lastgamma);
                ArrayList<Double> threshold = get_threshold_new(Constants.slide, Constants.slide, lastalpha, lastbeta, lastgamma);
                Collections.sort(threshold);
//                double threshold = get_dynamic_threshold(lastR, radius, Constants.W, Constants.slide,
//                        lastgamma, all_core_points_map.get(sIdx).get(lastR).size(), all_core_points_map.get(sIdx).get(radius).size());
                System.out.println("Computing for radius = " + radius);
                if (threshold.size() == 3) {
                    System.out.println("Threshold = " + threshold.get(0) + ", " + threshold.get(1) + ", " + threshold.get(2));
                } else {
                    System.out.println("Threshold = " + threshold.get(0));
                }
                if ((threshold.size() == 3
                        && (all_core_points_map.get(sIdx).get(radius).size() > threshold.get(2)
                        || (all_core_points_map.get(sIdx).get(radius).size() > threshold.get(0)
                        && all_core_points_map.get(sIdx).get(radius).size() < threshold.get(1))))
                        || (threshold.size() == 1 && all_core_points_map.get(sIdx).get(radius).size() > threshold.get(0))) {
                    System.out.println("Not creating cores for radius = " + radius);
                    needCore = false;
                }
            }
            if (needCore) {
                for (CorePoint c : all_core_points_map.get(sIdx).get(radius)) {
                    if (!checked_core.contains(c)) {
                        for (int i = 0; i < Constants.slide; i++) {
                            checked[i] = false;
                        }
                        for (Bin b : c.all_bins) {
                            if (b.data.containsKey(sIdx)) {
                                for (C_Data d : b.data.get(sIdx)) {
                                    checked[d.arrivalTime - all_slides.get(sIdx).get(0).arrivalTime] = true;
                                }
                            }
                        }
                        ArrayList<Double> all_r_in_core = new ArrayList<>(c.supported_r);
                        Double prob_r = all_r_in_core.get(0);
                        for (Double r : all_r_in_core) {
                            if (r > prob_r) {
                                prob_r = r;
                            }
                        }

                        for (CorePoint c2 : all_core_points_map.get(sIdx).get(prob_r)) {

                            if (c != c2) {
                                double distance = DistanceFunction.euclideanDistance(c, c2);

                                if (distance <= prob_r * 3) {
                                    checked = probCoreWithList(c, c2.getDataInRange(0, prob_r, sIdx), sIdx,
                                            all_r_in_core,
                                            prob_r,
                                            checked,
                                            all_slides.get(sIdx).get(0).arrivalTime);
                                }
                            }
                        }
                        checked_core.add(c);
                    }
                }
                //compute average T2r2
                double T2r2 = 0;
                double Thalfr2 = 0;
                double Tr2 = 0;
                for (CorePoint c : all_core_points_map.get(sIdx).get(radius)) {
                    T2r2 += c.getDataInRange(0, 2 * radius, sIdx).size();
                    Thalfr2 += c.getDataInRange(0, radius / 2, sIdx).size();
                    Tr2 += c.getDataInRange(0, radius, sIdx).size();
                }
                lastgamma = T2r2 / all_core_points_map.get(sIdx).get(radius).size();
                lastR = radius;
                lastalpha = T2r2 / Thalfr2;
                lastbeta = T2r2 / Tr2;
            } else {
                r_to_remove.add(radius);
            }
        }
        for (CorePoint c : newCores) {

            for (Double r : r_to_remove) {
                c.supported_r.remove(r);
            }
            if (!c.supported_r.isEmpty()) {
                all_distince_core_point.add(c);
            }

        }
        ArrayList<CorePoint> to_remove = new ArrayList<>();
        for (CorePoint c : all_distince_core_point) {
            if (c.supported_r.isEmpty()) {
                to_remove.add(c);
            }
        }
        for (CorePoint c : to_remove) {
            all_distince_core_point.remove(c);
        }
    }

    private boolean[] probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx,
            List<Double> all_r,
            double max_r,
            boolean[] checked, int start_time) {

        //create bins if not exists
        c.creatBins(all_r);
        if (candidates != null) {
            for (C_Data d2 : candidates) {
                if (!checked[d2.arrivalTime - start_time]) {
                    double distance = DistanceFunction.euclideanDistance(d2, c);
                    if (distance <= max_r * 2) //put d2 to correct bin
                    {
                        c.putDataToBin(d2, distance, sIdx);

                        for (Double r : all_r) {
                            if (distance <= r) {
                                if (d2.closeCoreMaps_R.get(r) == null) {
                                    d2.closeCoreMaps_R.put(r, c);
                                }
                            }
                        }

                    }
                    checked[d2.arrivalTime - start_time] = true;
                }
            }
        }
        return checked;
    }

    class Bin {

        public double min_val;
        public double max_val;
        public HashMap<Integer, ArrayList<C_Data>> data;

        public Bin(double min, double max) {
            this.min_val = min;
            this.max_val = max;
            this.data = new HashMap<>();
        }
    }

    class CorePoint extends C_Data {

        public ArrayList<Bin> all_bins = new ArrayList<>();
//        public C_Data orig_data;

        public HashMap<Double, Integer> totalHalfRPoints = new HashMap<>();

        public HashSet<Double> supported_r = new HashSet<>();
        public Integer supported_slides = -1;

        public void discountHalfPoints(int expiredSlideIndex, Double r) {
            if (totalHalfRPoints.containsKey(r)) {

                int t = 0;
                for (Bin b : all_bins) {
                    if (b.max_val <= r / 2) {
                        if (b.data.containsKey(expiredSlideIndex)) {
                            t += b.data.get(expiredSlideIndex).size();
                        }
                    }
                }
                totalHalfRPoints.put(r, totalHalfRPoints.get(r) - t);
            }
        }

        public void updateHalfRPoints(int newestSlideIndex, Double r) {

            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= r / 2) {
                    if (b.data.containsKey(newestSlideIndex)) {
                        t += b.data.get(newestSlideIndex).size();
                    }
                }
            }

            totalHalfRPoints.put(r, totalHalfRPoints.get(r) + t);

        }

        public void computeTotalHalfRPoints(Double r) {

            if (!totalHalfRPoints.containsKey(r)) {
                int t = 0;
                for (Bin b : all_bins) {
                    if (b.max_val <= r / 2) {
                        for (Integer idx : b.data.keySet()) {

                            t += b.data.get(idx).size();

                        }
                    }
                }
                totalHalfRPoints.put(r, t);
            }

        }

        public ArrayList<C_Data> getDataInRange(double min, double max) {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min && b.max_val <= max) {
                    for (Integer sIdx : b.data.keySet()) {
                        result.addAll(b.data.get(sIdx));
                    }
                }
            }
            return result;
        }

        public void creatBins(List<Double> all_r) {

            //create bins if not exists
            if (all_bins == null || all_bins.isEmpty()) {
                ArrayList<Double> all_r_threshold = new ArrayList<>();
                all_r_threshold.add(0.0);
                for (Double r : all_r) {
                    all_r_threshold.add(r / 2);
                    all_r_threshold.add(r);
                    all_r_threshold.add(3 * r / 2);
                    all_r_threshold.add(2 * r);
                }
                all_bins = new ArrayList<>();
                Collections.sort(all_r_threshold);
                for (int i = 1; i < all_r_threshold.size(); i++) {
                    Bin b = new Bin(all_r_threshold.get(i - 1), all_r_threshold.get(i));
                    all_bins.add(b);
                }
            }
        }

        public void putDataToBin(C_Data d, double distance, int sIdx) {
            for (Bin b : all_bins) {
                if ((b.min_val <= distance) && b.max_val >= distance) {
                    if (!b.data.containsKey(sIdx)) {
                        b.data.put(sIdx, new ArrayList<>());
                    }
                    b.data.get(sIdx).add(d);
                }
            }
        }

        public ArrayList<C_Data> getDataInRange(double min, double max, int sIdx) {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min && b.max_val <= max) {
                    if (b.data.containsKey(sIdx)) {
                        result.addAll(b.data.get(sIdx));
                    }
                }
            }
            return result;
        }

//        public HashMap<OD_Query, Integer> totalHalfRPoints = new HashMap<>();
        public int getTotalHalfRPoints(OD_Query q) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= q.R / 2) {
                    for (Integer idx : b.data.keySet()) {

                        t += b.data.get(idx).size();

                    }
                }
            }
            return t;
        }

        public int getTotalHalfRPoints(OD_Query q, int sIdx) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= q.R / 2) {
                    if (b.data.containsKey(sIdx)) {
                        t += b.data.get(sIdx).size();
                    }

                }
            }
            return t;
        }

        public CorePoint(C_Data d) {
            this.values = d.values;
            this.hashCode = d.hashCode;
        }

        private void reset() {
            for (Bin b : all_bins) {
                b.data.clear();
            }
            all_bins.clear();
            totalHalfRPoints.clear();
//            orig_data = null;
        }

    }

    class C_Data extends Data {

        public HashMap<Double, Integer> numSucceedingNeighbor = new HashMap<>();
//        public HashMap<Double, Integer> neighborCount = new HashMap<>();

        public HashMap<Integer, HashMap<Double, Integer>> pred_neighbor_count = new HashMap<>();

        public HashMap<Double, CorePoint> closeCoreMaps_halfR = new HashMap<>();
        public HashMap<Double, CorePoint> closeCoreMaps_R = new HashMap<>();

        public int lastProbRight = -1;
        public int lastProbLeft = -1;
        public int sIndex = -1;
        public ArrayList<OD_Query> safeInlierQueries = new ArrayList<>();

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;
            for (OD_Query q : all_queries) {
//                if (!neighborCount.containsKey(q.R)) {
//                    neighborCount.put(q.R, 0);
//                }
                if (!numSucceedingNeighbor.containsKey(q.R)) {
                    numSucceedingNeighbor.put(q.R, 0);
                }
            }
        }

        public C_Data() {

        }

        public HashMap<Double, Integer> neighborCounts() {
            HashMap<Double, Integer> neighborCount = new HashMap<>();
            for (Double r : all_r) {
                int count = get_neighbor_count(numSucceedingNeighbor, r);
                for (Integer sIdx : all_slides.keySet()) {
                    if (pred_neighbor_count.containsKey(sIdx)) {
                        count += get_neighbor_count(pred_neighbor_count.get(sIdx), r);
                    }
                }
                neighborCount.put(r, count);

            }
            return neighborCount;
        }

        public boolean isSupported() {
            return (closeCoreMaps_halfR != null)
                    || (!closeCoreMaps_R.isEmpty());
        }

        public ArrayList<OD_Query> getNeedProbQueries() {
            ArrayList<OD_Query> result = new ArrayList<>();
            HashMap<Double, Integer> neighborCount = neighborCounts();
            for (OD_Query q : all_queries) {
                if (!safeInlierQueries.contains(q) && neighborCount.get(q.R) < q.k) {

                    boolean need = true;
                    //check if connect to core point has more than k+1 neighbors in range R/2
                    if (closeCoreMaps_halfR.get(q.R) != null
                            && closeCoreMaps_halfR.get(q.R).totalHalfRPoints.get(q.R) != null
                            && closeCoreMaps_halfR.get(q.R).totalHalfRPoints.get(q.R) > q.k) {
                        need = false;
                    }
                    if (need) {
                        result.add(q);
                    }
                }
            }
            return result;
        }

        private void reset() {
            numSucceedingNeighbor.clear();
            pred_neighbor_count.clear();
//            neighborCount.clear();
            closeCoreMaps_R.clear();
            closeCoreMaps_halfR.clear();
            safeInlierQueries.clear();
        }

        private void refineNeighborStore() {
            //remove neighbor count = 0
            for (Integer s : pred_neighbor_count.keySet()) {
                HashMap<Double, Integer> neighborCountMap = pred_neighbor_count.get(s);
                for (Double r : all_r) {
                    if (neighborCountMap.containsKey(r) && neighborCountMap.get(r) <= 0) {
                        neighborCountMap.remove(r);
                    }
                }
            }
            for (Double r : all_r) {
                if (numSucceedingNeighbor.containsKey(r) && numSucceedingNeighbor.get(r) <= 0) {
                    numSucceedingNeighbor.remove(r);
                }
            }

            //keep only one smallest r with count >= k_max
            Integer k_max = all_k.get(all_k.size() - 1);
            for (int i = all_r.size() - 1; i >= 1; i--) {
                Double r = all_r.get(i);
                Double prev_r = all_r.get(i - 1);
                if (numSucceedingNeighbor.containsKey(r)
                        && numSucceedingNeighbor.containsKey(prev_r)
                        && numSucceedingNeighbor.get(r) >= k_max
                        && numSucceedingNeighbor.get(prev_r) >= k_max) {
                    numSucceedingNeighbor.remove(r);
                }

                for (Integer s : pred_neighbor_count.keySet()) {
                    HashMap<Double, Integer> neighborCountMap = pred_neighbor_count.get(s);
                    if (neighborCountMap.containsKey(r)
                            && neighborCountMap.containsKey(prev_r)
                            && neighborCountMap.get(r) >= k_max
                            && neighborCountMap.get(prev_r) >= k_max) {
                        neighborCountMap.remove(r);
                    }
                }
            }

            //remove duplicate neighbor count 
            for (int i = all_r.size() - 1; i >= 1; i--) {
                Double r = all_r.get(i);
                Double prev_r = all_r.get(i - 1);
                if (numSucceedingNeighbor.containsKey(r)
                        && numSucceedingNeighbor.containsKey(prev_r)
                        && Objects.equals(numSucceedingNeighbor.get(r), numSucceedingNeighbor.get(prev_r))) {
                    numSucceedingNeighbor.remove(r);
                }

                for (Integer s : pred_neighbor_count.keySet()) {
                    HashMap<Double, Integer> neighborCountMap = pred_neighbor_count.get(s);
                    if (neighborCountMap.containsKey(r)
                            && neighborCountMap.containsKey(prev_r)
                            && Objects.equals(neighborCountMap.get(r), neighborCountMap.get(prev_r))) {
                        neighborCountMap.remove(r);
                    }
                }
            }

        }

    }
}
