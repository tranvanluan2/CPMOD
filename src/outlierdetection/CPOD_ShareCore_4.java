/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.TreeMap;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author Luan Tran
 */
public class CPOD_ShareCore_4 {

    public static int currentTime;

    public static int expiredSlideIndex = -1;
    public static int newestSlide = -1;
    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();
    public HashSet<OD_Query> all_query_set = new HashSet<>();

    public ArrayList<Double> all_sorted_r;
    public ArrayList<Integer> all_sorted_k;
    public HashMap<Double, ArrayList<Integer>> all_r_w = new HashMap<>();
    public static HashMap<Double, ArrayList<OD_Query>> r_query_map = new HashMap<>();

    public static TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> rkw_map = new TreeMap<>();

    public static HashMap<Integer, HashMap<Double, ArrayList<CorePoint>>> all_core_points_map = new HashMap<>();

    public static ArrayList<CorePoint> all_distinct_core_point = new ArrayList<>();

//    public static HashMap<Double, ArrayList<Integer>> corepointstats = new HashMap<>();
    public HashMap<OD_Query, ArrayList<Data>> slide_process(ArrayList<Data> data, int _currentTime) {
        HashMap<OD_Query, ArrayList<Data>> result = new HashMap<>();
        for (OD_Query q : all_queries) {
            result.put(q, new ArrayList<>());
        }
        currentTime = _currentTime;
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
            ArrayList<C_Data> idx = all_slides.get(d.sIndex);
            if (idx != null) {
                idx.add(d);
            } else {
                idx = new ArrayList<>(Constants.slide);
                idx.add(d);
                all_slides.put(d.sIndex, idx);
            }
        }

        expiredSlideIndex = (currentTime - 1) / Constants.slide - Constants.W / Constants.slide;
        processExpiredData(expiredSlideIndex);
        long start = Utils.getCPUTime();
        for (int sIdx : slide_to_process) {
            selectAllCore(sIdx, all_sorted_r);
        }

        newestSlide = (currentTime - 1) / Constants.slide;

        if (currentTime == Constants.W) {
            for (CorePoint c : all_distinct_core_point) {
                for (Double r : all_sorted_r) {
                    if (r <= c.max_support_r) {
                        c.computeTotalHalfRPoints(r);
                    } else {
                        break;
                    }
                }
            }
        } else {
            for (Double r : all_sorted_r) {
                for (CorePoint c : all_core_points_map.get(newestSlide).get(r)) {
                    c.computeTotalHalfRPoints(r);
//                    if (c.totalHalfRPoints.get(r) != null) {
//                        c.updateHalfRPoints(r);
//                    } else {
//                        c.computeTotalHalfRPoints(r);
//                    }

                }
            }
        }
        System.out.println("Time for creating core points = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        start = Utils.getCPUTime();
        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                ArrayList<OD_Query> need_prob_queries = new ArrayList<>();
                for (OD_Query q : all_queries) {
                    //check with core points
                    if ((!d.closeCoreMaps_halfR.containsKey(q.R)

                            || !d.closeCoreMaps_halfR.get(q.R).enoughPointsInHalfR(q.R, q.W, q.k) //                            d.closeCoreMaps_halfR.get(q.R).totalHalfRPoints.get(q.R).get(q.W) <= q.k
                            )
                            && !d.safeInlierQueries.contains(q)) {
                        if (currentTime - q.W < d.arrivalTime) {
                            need_prob_queries.add(q);
                        }
                    }
                }

                //compute k max and r-w map
                Integer k_max = all_sorted_k.get(all_sorted_k.size() - 1);
                TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> r_w_map = new TreeMap<>();
                for (OD_Query q : need_prob_queries) {
                    if (!r_w_map.containsKey(q.R)) {
                        r_w_map.put(q.R, new TreeMap<>());
                    }
                    if (!r_w_map.get(q.R).containsKey(q.W)) {
                        r_w_map.get(q.R).put(q.W, new ArrayList<>());
                    }
                    r_w_map.get(q.R).get(q.W).add(q);

                }
                for (Double r : r_w_map.keySet()) {
                    ArrayList<Integer> all_w = new ArrayList<>(r_w_map.get(r).keySet());
                    int total_neighbor = 0;
                    int total_suc_neighbor = 0;
                    for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {

                        if (d.lastProbRadius.containsKey(idx) && d.lastProbRadius.get(idx) >= r) {
                            int n = d.get_neighbor_count(d.neighborCount.get(idx), r);
                            total_neighbor += n;
                            if (idx >= d.sIndex) {
                                total_suc_neighbor += n;
                            }
                        } else {
                            int[] ns = probe_slide(d, idx, r, k_max, total_neighbor, total_suc_neighbor);
                            total_neighbor = ns[0];
                            total_suc_neighbor = ns[1];
                        }

                        if (total_neighbor >= k_max) {
                            if (total_suc_neighbor >= k_max) {
                                for (Integer w : all_w) {
                                    if (idx >= newestSlide - all_w.get(0) / Constants.slide + 1) {
                                        d.safeInlierQueries.addAll(r_w_map.get(r).get(w));
                                    }
                                }
                            }
                            break;
                        } else {
                            if (idx == newestSlide - all_w.get(0) / Constants.slide + 1) {
                                for (OD_Query q : r_w_map.get(r).get(all_w.get(0))) {
                                    if (total_neighbor < q.k) {
                                        if (!result.containsKey(q)) {
                                            result.put(q, new ArrayList<>());
                                        }
                                        result.get(q).add(d);
                                    }
                                }
                                all_w.remove(0);
                                if (all_w.isEmpty()) {
                                    break;
                                }
                            }
                        }

                    }

                }
            }
        }

        double timeForFindNeighbors = (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        System.out.println("Time for finding neighbor = " + timeForFindNeighbors);
        return result;
    }

    private int[] probe_slide(C_Data d, int sIdx, Double r, int k_max, int cur_neighbor, int cur_sucneighbor) {
        if (!d.neighborCount.containsKey(sIdx)) {
            d.neighborCount.put(sIdx, new TreeMap<>());
        }
        //scan possible points 
        int min_arrival_time = all_slides.get(sIdx).get(0).arrivalTime;
        boolean[] checked = null;
        int case_ = 0;
        ArrayList<Bin> possibleCandidates = new ArrayList<>();
        Double lastR = -1.0;
        if (!d.lastProbRadius.containsKey(sIdx) || d.lastProbCore.get(sIdx) == null) {
            //find close core
            ResultFindCore rf = findCloseCore(d, sIdx, r);
            double distance = rf.getDistance();
            ArrayList<CorePoint> cores = rf.getCore();
            if (cores != null) {
                if (distance <= r / 2) {
                    CorePoint c = cores.get(0);
                    //grab close neighbor in range R/2 of c
                    possibleCandidates.addAll(c.getDataInRange(0, 1.5 * r, sIdx));
                    d.lastProbCore.put(sIdx, c);
                } else if (distance <= r) {
                    possibleCandidates.addAll(cores.get(0).getDataInRange(0, 2 * r, sIdx));
                    d.lastProbCore.put(sIdx, cores.get(0));
                } else if (distance <= r * 2) {
                    case_ = 1;
                    checked = new boolean[Constants.slide];
                    for (CorePoint c : cores) {
                        possibleCandidates.addAll(c.getDataInRange(0, r, sIdx));
                    }
                    d.lastProbCore.put(sIdx, null);
                }

            }
        } else if (d.lastProbRadius.get(sIdx) < r && d.lastProbCore.get(sIdx) != null) {

            CorePoint c = d.lastProbCore.get(sIdx);
            double distance = DistanceFunction.euclideanDistance(d, c);
            lastR = d.lastProbRadius.get(sIdx);
            double lowerbound = lastR / 2;
            if (distance <= lastR / 2) {
                lowerbound = lastR / 2;
            }
            //else if (distance <= lastR) {
            //  lowerbound = 2 * lastR;
            //}
            if (distance <= r / 2) {
                //grab close neighbor in range R/2 of c
                possibleCandidates.addAll(c.getDataInRange(lowerbound, 1.5 * r, sIdx));
            } else if (distance <= r) {
                possibleCandidates.addAll(c.getDataInRange(lowerbound, 2 * r, sIdx));
            }

        }

        outer:
        for (Bin ps : possibleCandidates) {
            for (int t = 0; t < ps.data.get(sIdx).size(); t++) {
                C_Data d2 = ps.data.get(sIdx).get(t);
                if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {
                    double dist = DistanceFunction.euclideanDistance(d2, d);
                    if (dist > lastR && dist <= r) {
                        Double matchR = getMatchRadius(dist, all_sorted_r);
                        if (!d.neighborCount.get(sIdx).containsKey(matchR)) {
                            d.neighborCount.get(sIdx).put(matchR, 1);
                        } else {
                            d.neighborCount.get(sIdx).put(matchR, d.neighborCount.get(sIdx).get(matchR) + 1);
                        }
                        cur_neighbor += 1;
                        if (d.sIndex <= sIdx) {
                            cur_sucneighbor += 1;
                        }
                        //check condition
                        if ((cur_neighbor >= k_max && d.sIndex > sIdx) || (d.sIndex <= sIdx && cur_sucneighbor >= k_max)) {
                            break outer;
                        }
                    }

                    if (case_ == 1) {
                        checked[d2.arrivalTime - min_arrival_time] = true;
                    }

                }
            }
        }

        d.lastProbRadius.put(sIdx, r);
        return new int[]{cur_neighbor, cur_sucneighbor};

    }

//    private ArrayList<OD_Query> probe_r_k(C_Data d, int newestSlide, int oldestSlide,
//            ArrayList<OD_Query> need_prob_queries) {
//        ArrayList<OD_Query> outlierQueries = new ArrayList<>();
//        TreeMap<Double, ArrayList<OD_Query>> queryMap = convertToMap(need_prob_queries);
//        ArrayList<Double> all_r = new ArrayList<>(queryMap.keySet());
//        int k_max = -1;
//        for (OD_Query q : need_prob_queries) {
//            if (q.k > k_max) {
//                k_max = q.k;
//            }
//        }
//        while (!all_r.isEmpty()) {
//            int selected_idx = 0;
//            Double r = all_r.get(selected_idx);
//            int k = queryMap.get(r).get(0).k; //max k
//            boolean isInlier = false;
//            boolean safeInlier = false;
//            int cur_neighbor = d.neighborCount(r, oldestSlide);
//            int cur_suc_neighbor = d.sucNeighborCount(r);
//            if (cur_neighbor >= k) {
//                isInlier = true;
//                if (cur_suc_neighbor >= k) {
//                    safeInlier = true;
//                }
//            }
//
//            //probe right 
//            if (!isInlier) {
//                for (int sIdx = d.sIndex; sIdx <= newestSlide; sIdx++) {
//                    if (d.lastProbRadius.containsKey(sIdx) && d.lastProbRadius.get(sIdx) >= r) {
//                        continue;
//                    }
//                    int[] ns = probe_slide(d, sIdx, r, k_max, cur_neighbor, cur_suc_neighbor);
//                    cur_neighbor = ns[0];
//                    cur_suc_neighbor = ns[1];
//                    //check inlier 
//                    if (cur_neighbor >= k) {
//                        isInlier = true;
//                        //check safe inlier
//                        if (cur_suc_neighbor >= k) {
//                            safeInlier = true;
//                            for (OD_Query q : queryMap.get(r)) {
//                                d.safeInlierQueries.add(q);
//                            }
//                        }
//                        break;
//                    }
//
//                }
//            }
//            //probe left if not inlier
//            if (!isInlier) {
//                for (int sIdx = d.sIndex - 1; sIdx >= oldestSlide; sIdx--) {
//                    if (d.lastProbRadius.containsKey(sIdx) && d.lastProbRadius.get(sIdx) >= r) {
//                        continue;
//                    }
//                    int[] ns = probe_slide(d, sIdx, r, k_max, cur_neighbor, cur_suc_neighbor);
//                    cur_neighbor = ns[0];
//                    cur_suc_neighbor = ns[1];
//                    //check inlier 
//                    if (cur_neighbor >= k) {
//                        isInlier = true;
//                        break;
//                    }
//
//                }
//            }
//
//            if (isInlier) {
//                if (safeInlier) {
//
//                    for (Double largerR : all_r.subList(selected_idx + 1, all_r.size())) {
//
//                        for (OD_Query q : queryMap.get(largerR)) {
//                            if (q.k <= k) {
//                                d.safeInlierQueries.add(q);
//                            }
//                        }
//
//                    }
//                }
//                for (int r_idx = all_r.size() - 1; r_idx > selected_idx; r_idx--) {
//                    Double largerR = all_r.get(r_idx);
//                    if (k >= queryMap.get(largerR).get(0).k) {
//                        all_r.remove(r_idx);
//                    }
//                }
//
//            } else {
//                for (OD_Query q : queryMap.get(r)) {
//                    if (q.k >= cur_neighbor) {
//                        outlierQueries.add(q);
//                    } else {
//                        break;
//                    }
//                }
//                for (int r_idx = selected_idx - 1; r_idx >= 0; r_idx--) {
//                    Double smallerR = all_r.get(r_idx);
//                    int smallestK = queryMap.get(smallerR).get(queryMap.get(smallerR).size() - 1).k;
//                    if (smallerR < r && smallestK >= k) {
//
//                        for (OD_Query q : queryMap.get(smallerR)) {
//
//                            outlierQueries.add(q);
//
//                        }
//                        all_r.remove(r_idx);
//                    }
//                }
//            }
//            all_r.remove(r);
//
//        }
//        return outlierQueries;
//    }
//    private ArrayList<OD_Query> probe(C_Data d, int newestSlide, ArrayList<OD_Query> need_prob_queries) {
//        ArrayList<OD_Query> outlierQueries = new ArrayList<>();
//        TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> queryMapRKW
//                = convertToRKW(need_prob_queries);
//
//        return outlierQueries;
//    }
    public Double getMatchRadius(double distance, ArrayList<Double> all_rs) {
        if (distance <= all_rs.get(0)) {
            return all_rs.get(0);
        }
        if (distance > all_rs.get(all_rs.size() - 1)) {
            return null;
        }
        int low = 0;
        int high = all_rs.size() - 1;
        while (low < high - 1) {
            int mid = low + (high - low) / 2;
            if (all_rs.get(mid) < distance) {
                low = mid;
            } else {
                high = mid;
            }
        }
        return all_rs.get(high);
    }

//    private TreeMap<Double, ArrayList<OD_Query>> convertToMap(ArrayList<OD_Query> need_prob_queries) {
//        TreeMap<Double, ArrayList<OD_Query>> result = new TreeMap<>();
//        for (OD_Query q : need_prob_queries) {
//            if (!result.containsKey(q.R)) {
//                result.put(q.R, new ArrayList<>());
//            }
//            result.get(q.R).add(q);
//        }
//        return result;
//    }
    private TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> convertToRKW(ArrayList<OD_Query> need_prob_queries) {
        TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> results = new TreeMap<>();
        for (OD_Query q : need_prob_queries) {
            if (!results.containsKey(q.R)) {
                results.put(q.R, new TreeMap<>());
            }
            if (!results.get(q.R).containsKey(q.k)) {
                results.get(q.R).put(q.k, new TreeMap<>());
            }
            if (!results.get(q.R).get(q.k).containsKey(q.W)) {
                results.get(q.R).get(q.k).put(q.W, new ArrayList<>());
            }
            results.get(q.R).get(q.k).get(q.W).add(q);
        }
        return results;

    }

//    private TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> convertToRKW(HashSet<OD_Query> need_prob_queries) {
//        TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> results = new TreeMap<>();
//        for (OD_Query q : need_prob_queries) {
//            if (!results.containsKey(q.R)) {
//                results.put(q.R, new TreeMap<>());
//            }
//            if (!results.get(q.R).containsKey(q.k)) {
//                results.get(q.R).put(q.k, new TreeMap<>());
//            }
//            if (!results.get(q.R).get(q.k).containsKey(q.W)) {
//                results.get(q.R).get(q.k).put(q.W, new ArrayList<>());
//            }
//            results.get(q.R).get(q.k).get(q.W).add(q);
//        }
//        return results;
//
//    }
//    private TreeMap<Integer, ArrayList<OD_Query>> convertToMapR_W(ArrayList<OD_Query> need_prob_queries) {
//        TreeMap<Integer, ArrayList<OD_Query>> result = new TreeMap<>();
//        for (OD_Query q : need_prob_queries) {
//            if (!result.containsKey(q.W)) {
//                result.put(q.W, new ArrayList<>());
//            }
//            result.get(q.W).add(q);
//        }
////        for (Integer w : result.keySet()) {
////            Collections.sort(result.get(w), new RKComparator());
////        }
//        return result;
//    }
    private ResultFindCore findCloseCore(C_Data d, int slideIndex, double radius) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR.get(radius) != null
                && d.closeCoreMaps_halfR.get(radius).supported_slides.contains(slideIndex)) {
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

    public void add_query(OD_Query q) {
        all_queries.add(q);
        all_query_set.add(q);
        //sort queries in increasing R order and reverse k
        Collections.sort(all_queries, new RKComparator());
        if (!r_query_map.containsKey(q.R)) {
            r_query_map.put(q.R, new ArrayList<>());
        }
        r_query_map.get(q.R).add(q);
        Collections.sort(r_query_map.get(q.R), new RKComparator());
        get_unique_r();
        get_unique_k();
        rkw_map = convertToRKW(all_queries);

        if (!all_r_w.containsKey(q.R)) {
            all_r_w.put(q.R, new ArrayList<>());
        }
        if (!all_r_w.get(q.R).contains(q.W)) {
            all_r_w.get(q.R).add(q.W);
        }
        Collections.sort(all_r_w.get(q.R));
    }

    public void get_unique_k() {
        all_sorted_k = new ArrayList<>();
        for (OD_Query q : all_queries) {
            if (!all_sorted_k.contains(q.k)) {
                all_sorted_k.add(q.k);
            }
        }
        Collections.sort(all_sorted_k);
    }

    public void get_unique_r() {
        all_sorted_r = new ArrayList<>();
        for (OD_Query q : all_queries) {
            if (!all_sorted_r.contains(q.R)) {
                all_sorted_r.add(q.R);
            }
        }
        Collections.sort(all_sorted_r);
    }

    private void processExpiredData(int expiredSlideIndex) {
        all_slides.remove(expiredSlideIndex);
        for (int c_idx = all_distinct_core_point.size() - 1; c_idx >= 0; c_idx--) {
            CorePoint c = all_distinct_core_point.get(c_idx);
            for (Double r : all_sorted_r) {
                if (r <= c.max_support_r) {
                    c.discountHalfPoints(r);
                } else {
                    break;
                }
            }
            for (Bin b : c.all_bins) {
//                if (b.data.containsKey(expiredSlideIndex)) {
                b.data.remove(expiredSlideIndex);
//                }
            }
            c.supported_slides.remove(expiredSlideIndex);
            if (c.supported_slides.isEmpty()) {
                all_distinct_core_point.remove(c_idx);
            }
        }
        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
//                if (d.neighborCount.containsKey(expiredSlideIndex)) {
                d.neighborCount.remove(expiredSlideIndex);
//                }
            }
        }
//        if (all_core_points_map.containsKey(expiredSlideIndex)) {
//            all_core_points_map.get(expiredSlideIndex).clear();
//        }
        all_core_points_map.remove(expiredSlideIndex);
    }

    private void selectAllCore(int sIdx, ArrayList<Double> all_sorted_r) {
        HashSet<CorePoint> indexed_cores = new HashSet<>();
        ArrayList<C_Data> data = all_slides.get(sIdx);
        ArrayList<CorePoint> current_cores = new ArrayList<>();
        for (int r_idx = all_sorted_r.size() - 1; r_idx >= 0; r_idx--) {
            Double r = all_sorted_r.get(r_idx);
            for (C_Data d : data) {
                if (d.closeCoreMaps_halfR.get(r) == null && d.closeCoreMaps_R.get(r) == null) {
                    //check with current cores
                    for (CorePoint c : current_cores) {
                        double dist = DistanceFunction.euclideanDistance(d, c);
                        if (dist <= r) {
                            //connect d to c
                            d.connectCore(c, dist, all_sorted_r.subList(0, r_idx + 1));
                            if (!indexed_cores.contains(c)) {
                                c.putDataToBin(d, dist, sIdx);
                            }
                            break;
                        }
                    }
                }
                //if d is not linked to a core, find in history
                if (d.closeCoreMaps_halfR.get(r) == null && d.closeCoreMaps_R.get(r) == null) {
                    for (CorePoint c : all_distinct_core_point) {
                        if (c.max_support_r >= r) {
                            double dist = DistanceFunction.euclideanDistance(d, c);
                            if (dist <= r) {
                                d.connectCore(c, dist, all_sorted_r.subList(0, r_idx + 1));
                                current_cores.add(c);
                                if (!indexed_cores.contains(c)) {
                                    c.putDataToBin(d, dist, sIdx);
                                }
                                break;
                            }
                        }
                    }

                }

                //if d is not linked to any core in history
                if (d.closeCoreMaps_halfR.get(r) == null && d.closeCoreMaps_R.get(r) == null) {
                    //creat a core from d
                    CorePoint c = new CorePoint(d);
                    c.creatBins(all_sorted_r.subList(0, r_idx + 1));
                    c.max_support_r = r;
                    current_cores.add(c);
                    c.putDataToBin(d, 0, sIdx);
                    d.connectCore(c, 0, all_sorted_r.subList(0, r_idx + 1));
                    all_distinct_core_point.add(c);
                }
            }

            //index data points around cores
            for (CorePoint c : current_cores) {
                if (!indexed_cores.contains(c)) {
                    boolean[] checked = new boolean[Constants.slide];
                    for (Bin b : c.all_bins) {
                        if (b.data.containsKey(sIdx)) {
                            for (C_Data d : b.data.get(sIdx)) {
                                checked[d.arrivalTime - all_slides.get(sIdx).get(0).arrivalTime] = true;
                            }
                        }
                    }
                    for (CorePoint c2 : current_cores) {
                        if (c != c2) {
                            double distance = DistanceFunction.euclideanDistance(c, c2);
                            if (distance <= r * 3) {
                                for (Bin candidate : c2.getDataInRange(0, r, sIdx)) {
                                    checked = probCoreWithList(c, candidate.data.get(sIdx), sIdx,
                                            all_sorted_r.subList(0, r_idx + 1),
                                            checked,
                                            all_slides.get(sIdx).get(0).arrivalTime);
                                }
                            }
                        }
                    }
                    indexed_cores.add(c);
                }
            }

            //add to all_corepoint map
            if (!all_core_points_map.containsKey(sIdx)) {
                all_core_points_map.put(sIdx, new HashMap<>());
            }
            all_core_points_map.get(sIdx).put(r, (ArrayList<CorePoint>) current_cores.clone());
            for (CorePoint c : current_cores) {
                c.supported_slides.add(sIdx);
            }
        }

    }

    private boolean[] probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx,
            List<Double> all_r,
            boolean[] checked, int start_time) {
        Double max_r = all_r.get(all_r.size() - 1);
        if (candidates != null) {
            for (C_Data d2 : candidates) {
                if (!checked[d2.arrivalTime - start_time]) {
                    double distance = DistanceFunction.euclideanDistance(d2, c);
                    if (distance <= max_r * 2) //put d2 to correct bin
                    {
                        c.putDataToBin(d2, distance, sIdx);
                        d2.connectCore(c, distance, all_r);
                    }
                    checked[d2.arrivalTime - start_time] = true;
                }
            }
        }
        return checked;
    }

    class RWComparator implements Comparator<OD_Query> {

        @Override
        public int compare(OD_Query o1, OD_Query o2) {
            if (o1.R > o2.R) {
                return 2;
            } else if (o1.R == o2.R) {
                if (o1.W < o2.W) {
                    return 1;
                } else if (o1.W > o2.W) {
                    return -1;
                } else {
                    return 0;
                }
            } else {
                return -2;
            }
        }

    }

    class RKComparator implements Comparator<OD_Query> {

        @Override
        public int compare(OD_Query o1, OD_Query o2) {
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
        }

    }

    class CorePoint extends C_Data {

        public ArrayList<Bin> all_bins = new ArrayList<>();
//        public HashMap<Double, Integer> totalHalfRPoints = new HashMap<>();

        public HashMap<Double, TreeMap<Integer, Integer>> totalHalfRPoints = new HashMap<>();
        public Double max_support_r;
        public HashSet<Integer> supported_slides = new HashSet<>();
//        public HashMap<Double, Integer> lastHalfR = new HashMap<>();

        public boolean enoughPointsInHalfR(Double r, Integer w, int k_max) {
            if (totalHalfRPoints.get(r).get(w) != null) {
                return totalHalfRPoints.get(r).get(w) > k_max;
            } else {
                return true;
            }
//            {
//                for(Integer w2: totalHalfRPoints.get(r).descendingKeySet()){
//                    if(w2 <= w){
//                        return totalHalfRPoints.get(r).get(w2) > k_max;
//                    }
//                }
//            }

//            return false;
        }

        public void discountHalfPoints(Double r) {
            if (totalHalfRPoints.containsKey(r)) {

                for (Integer w : totalHalfRPoints.get(r).keySet()) {
                    int expiredSlide = newestSlide - w / Constants.slide;
                    int t = 0;
                    for (Bin b : all_bins) {
                        if (b.max_val <= r / 2) {
                            if (b.data.containsKey(expiredSlide)) {
                                t += b.data.get(expiredSlide).size();
                            }
                        } else {
                            break;
                        }
                    }

                    totalHalfRPoints.get(r).put(w, totalHalfRPoints.get(r).get(w) - t);
                }

            }
        }

        public void updateHalfRPoints(Double r) {

            if (totalHalfRPoints.containsKey(r)) {
                int t = 0;
                for (Bin b : all_bins) {
                    if (b.max_val <= r / 2) {
                        if (b.data.containsKey(newestSlide)) {
                            t += b.data.get(newestSlide).size();
                        }
                    } else {
                        break;
                    }
                }
                for (Integer w : totalHalfRPoints.get(r).keySet()) {
                    totalHalfRPoints.get(r).put(w, totalHalfRPoints.get(r).get(w) + t);
                }
            }
        }

        public void computeTotalHalfRPoints(Double r) {
            ArrayList<Integer> all_w = new ArrayList<>(all_r_w.get(r));
            TreeMap<Integer, Integer> numHalfR = new TreeMap<>();
            int t = 0;
            for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                for (Bin b : all_bins) {
                    if (b.max_val <= r / 2) {
                        if (b.data.containsKey(idx)) {
                            t += b.data.get(idx).size();
                        }
                    } else {
                        break;
                    }
                }
//                System.out.println("idx = "+idx);
//                System.out.println("other slide = "+(newestSlide - all_w.get(0) / Constants.slide + 1));
                while (!all_w.isEmpty() && idx == newestSlide - all_w.get(0) / Constants.slide + 1) {

                    numHalfR.put(all_w.get(0), t);

                    all_w.remove(0);

                    //break if t is larger than k_max
                }
                if (t > all_sorted_k.get(all_sorted_k.size() - 1)) {
                    break;
                }
            }

            totalHalfRPoints.put(r, numHalfR);
        }

//        public int computeTotalHalfRPoints(Double r, int idx) {
//            int t = 0;
//            for (Bin b : all_bins) {
//                if (b.max_val <= r / 2) {
//                    if (b.data.get(idx) != null) {
//                        t += b.data.get(idx).size();
//                    }
//
//                } else {
//                    break;
//                }
//            }
//            return t;
//        }
        public ArrayList<Bin> getDataInRange(double min, double max) {
            ArrayList<Bin> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min && b.max_val <= max) {
                    result.add(b);
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

        public ArrayList<Bin> getDataInRange(double min, double max, int sIdx) {
            ArrayList<Bin> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min) {
                    if (b.max_val <= max) {
                        if (b.data.containsKey(sIdx)) {
                            result.add(b);
                        }
                    }
                    if (b.min_val > max) {
                        break;
                    }
                }

            }
            return result;
        }

        public int getTotalHalfRPoints(Double r) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= r / 2) {
                    for (Integer idx : b.data.keySet()) {

                        t += b.data.get(idx).size();

                    }
                }
            }
            return t;
        }

        public int getTotalHalfRPoints(Double r, int sIdx) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= r / 2) {
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

    class C_Data extends Data {

        public HashMap<Integer, TreeMap<Double, Integer>> neighborCount = new HashMap<>();
        public HashMap<Double, CorePoint> closeCoreMaps_halfR = new HashMap<>();
        public HashMap<Double, CorePoint> closeCoreMaps_R = new HashMap<>();
        public HashMap<Integer, Double> lastProbRadius = new HashMap<>();
        public HashMap<Integer, CorePoint> lastProbCore = new HashMap<>();

        public int sIndex = -1;
        public HashSet<OD_Query> safeInlierQueries = new HashSet<>();

        public void connectCore(CorePoint c, double distance, List<Double> all_r) {
            for (int r_idx = all_r.size() - 1; r_idx >= 0; r_idx--) {
                Double r = all_r.get(r_idx);
                if (distance <= r) {
                    if (distance <= r / 2) {
                        closeCoreMaps_halfR.put(r, c);
                    } else {
                        closeCoreMaps_R.put(r, c);
                    }

                } else {
                    break;
                }
            }
        }

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.sIndex = (arrivalTime - 1) / Constants.slide;
        }

        public C_Data() {

        }

        private Integer get_neighbor_count(TreeMap<Double, Integer> neighbor_count_map, Double r) {
            int count = 0;
            for (Double ri : neighbor_count_map.keySet()) {
                if (ri <= r) {
                    count += neighbor_count_map.get(ri);
                } else {
                    break;
                }
            }

            return count;
        }

        private HashMap<Double, Integer> get_neighbor_count(TreeMap<Double, Integer> neighbor_count_map, ArrayList<Double> r_list) {
            HashMap<Double, Integer> result = new HashMap<>();
            int count = 0;
            int cur_idx = 0;
            for (Double ri : neighbor_count_map.keySet()) {
                count += neighbor_count_map.get(ri);
                if (Objects.equals(ri, r_list.get(cur_idx))) {
                    result.put(r_list.get(cur_idx), count);
                    cur_idx += 1;
                }

            }
            return result;
        }

        public HashMap<Double, Integer> neighborCounts(ArrayList<Double> r_list) {
            HashMap<Double, Integer> result = new HashMap<>();
            for (Integer sIdx : neighborCount.keySet()) {
                HashMap<Double, Integer> cur_neigh = get_neighbor_count(neighborCount.get(sIdx), r_list);
                for (Double r : r_list) {
                    if (!result.containsKey(r)) {
                        result.put(r, cur_neigh.get(r));
                    } else {
                        result.put(r, result.get(r) + cur_neigh.get(r));
                    }
                }
            }
            return result;

        }

        public int sucNeighborCount(Double r) {
            int count = 0;
            for (Integer sIdx : neighborCount.keySet()) {
                if (sIdx >= this.sIndex) {
                    count += get_neighbor_count(neighborCount.get(sIdx), r);
                }
            }
            return count;
        }

        public int neighborCount(Double r, int smallestSlideIndex) {
            int count = 0;
            for (Integer sIdx : neighborCount.keySet()) {
                if (sIdx >= smallestSlideIndex) {
                    count += get_neighbor_count(neighborCount.get(sIdx), r);
                }

            }
            return count;
        }

        public int neighborCount(Double r) {
            int count = 0;
            for (Integer sIdx : neighborCount.keySet()) {
                count += get_neighbor_count(neighborCount.get(sIdx), r);

            }
            return count;
        }

        public HashMap<Double, Integer> neighborCounts() {
            return neighborCounts(all_sorted_r);
        }

    }
}
