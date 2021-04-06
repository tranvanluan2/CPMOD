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
public class CPMOD_SD {

    public static int currentTime;

    public static int expiredSlideIndex = -1;
    public static int newestSlide = -1;
    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();
    public HashSet<OD_Query> all_query_set = new HashSet<>();

    public ArrayList<Double> all_sorted_r;
    public ArrayList<Integer> all_sorted_k;
    public HashMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> all_r_w = new HashMap<>();
    public static HashMap<Double, ArrayList<OD_Query>> r_query_map = new HashMap<>();

    public static TreeMap<Double, TreeMap<Integer, TreeMap<Integer, ArrayList<OD_Query>>>> rkw_map = new TreeMap<>();

    public static HashMap<Integer, HashMap<Double, ArrayList<CorePoint>>> all_core_points_map = new HashMap<>();

    public static ArrayList<CorePoint> all_distinct_core_point = new ArrayList<>();
    public static ArrayList<Integer> all_s = new ArrayList<>();

    public static int cur_slide_idx = -1;

    public static ArrayList<Integer> current_slides = new ArrayList<>();

    public static int last_report_slide = -1;
    public static boolean isOverlappingSlide = false;
    public static int last_report_time = -1;

    public HashMap<OD_Query, ArrayList<Data>> slide_process(ArrayList<Data> data, int _currentTime) {
        HashMap<OD_Query, ArrayList<Data>> result = new HashMap<>();
        for (OD_Query q : all_queries) {
            result.put(q, new ArrayList<>());
        }
        currentTime = _currentTime;
        boolean need_report = false;

        if (currentTime == Constants.W) {
            need_report = true;
            isOverlappingSlide = false;
            for (Data o : data) {
                if (isStartOfSlide(o.arrivalTime)) {
                    cur_slide_idx++;
                    newestSlide = cur_slide_idx;
                    current_slides.add(cur_slide_idx);
                }

                C_Data d = new C_Data(o);
                d.sIndex = cur_slide_idx;
                ArrayList<C_Data> idx = all_slides.get(d.sIndex);
                if (idx != null) {
                    idx.add(d);
                } else {
                    idx = new ArrayList<>(Constants.slide);
                    idx.add(d);
                    all_slides.put(d.sIndex, idx);
                }

            }

        } else {
            if (isStartOfSlide(data.get(0).arrivalTime)) {
                cur_slide_idx++;
                newestSlide = cur_slide_idx;
//                    if (current_slides.get(current_slides.size() - 1) < cur_slide_idx) {
                current_slides.add(cur_slide_idx);
//                    }
            } else {
                if (cur_slide_idx == last_report_slide) {
                    isOverlappingSlide = true;
                }
            }

            for (Data o : data) {

                C_Data d = new C_Data(o);
                d.sIndex = cur_slide_idx;
                ArrayList<C_Data> idx = all_slides.get(d.sIndex);
                if (idx != null) {
                    idx.add(d);
                } else {
                    idx = new ArrayList<>();
                    idx.add(d);
                    all_slides.put(d.sIndex, idx);
                }
            }

            Data o = data.get(data.size() - 1);
            if (isStartOfSlide(o.arrivalTime - Constants.W + 1)) {
                need_report = true;
            }

        }

//        for (Integer idx : all_slides.keySet()) {
//            System.out.println("Idx = " + idx + ", first arrival time = " + all_slides.get(idx).get(0).arrivalTime);
//        }
        if (need_report) {
            if (data.size() != Constants.W) {
                for (int i = 0; i < current_slides.size(); i++) {
                    if (all_slides.get(current_slides.get(i)).get(0).arrivalTime < currentTime + 1 - Constants.W) {
                        expiredSlideIndex = current_slides.get(i);
                    } else {
                        break;
                    }
                }
                while (current_slides.get(0) <= expiredSlideIndex) {
                    current_slides.remove(0);
                }

            }

//            ArrayList<C_Data> first_slide = all_slides.get(current_slides.get(0));
//            if (first_slide.get(0).arrivalTime != currentTime + 1 - Constants.W) {
//                System.out.println("Not right!");
//                System.out.println("First data arrival time = " + first_slide.get(0).arrivalTime);
//                System.out.println("Current time + 1- w = " + (currentTime + 1 - Constants.W));
//
//            }
            processExpiredData(expiredSlideIndex);
            long start = Utils.getCPUTime();
            if (data.size() == Constants.W) {
                for (int sIdx : current_slides) {
                    selectAllCore(sIdx, all_sorted_r);
                }
//                for (CorePoint c : all_distinct_core_point) {
//                    for (Double r : all_sorted_r) {
//                        if (r <= c.max_support_r) {
//                            c.computeTotalHalfRPoints(r);
//                        } else {
//                            break;
//                        }
//                    }
//                }
            } else {
                newestSlide = current_slides.get(current_slides.size() - 1);
                int start_slide = last_report_slide + 1;
                if (isOverlappingSlide) {
                    start_slide = last_report_slide;
                }
                for (int idx = start_slide; idx <= newestSlide; idx++) {
                    selectAllCore(idx, all_sorted_r);

                }
            }

            System.out.println("Time for creating core points = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

            ArrayList<OD_Query> reportingQueries = new ArrayList<>();
            for (OD_Query q : all_queries) {
                if ((currentTime - Constants.W) / q.S * q.S == (currentTime - Constants.W)) {
                    reportingQueries.add(q);
                }
            }

            HashMap<OD_Query, HashSet<CorePoint>> goodCores = new HashMap<>();
            for (OD_Query q : reportingQueries) {
                goodCores.put(q, new HashSet<>());
            }
            for (Integer sIdx : current_slides) {
                for (C_Data d : all_slides.get(sIdx)) {
                    //check valid and safe inlier queries
                    ArrayList<OD_Query> valid_queries = new ArrayList<>();
                    for (OD_Query q : reportingQueries) {
                        if (currentTime - q.W < d.arrivalTime && !d.safeInlierQueries.contains(q)) {
                            valid_queries.add(q);
                        }
                    }

                    ArrayList<OD_Query> need_prob_queries = new ArrayList<>();
                    CorePoint c = null;
                    HashSet<OD_Query> checked_queries = new HashSet<>();
                    for (int q_idx = 0; q_idx < valid_queries.size(); q_idx++) {
                        OD_Query q = valid_queries.get(q_idx);
                        if (checked_queries.contains(q)) {
                            continue;
                        }
                        //check core 
                        if (c == null) {
                            c = findCorePointInRangeR2(d, q.R, sIdx);
                        }
                        if (c != null && goodCores.get(q).contains(c)) {

                        } else if (c == null || !c.enoughPointsInHalfR(q.R, q.W, q.k)) {
                            need_prob_queries.add(q);
                        } else {
                            goodCores.get(q).add(c);
                            //mark check dominated queries
                            for (int q_idx2 = valid_queries.size() - 1; q_idx2 > q_idx; q_idx2--) {
                                OD_Query q2 = valid_queries.get(q_idx2);
                                if (q2.R >= q.R && q2.k <= q.k && q2.W >= q.W) {
                                    checked_queries.add(q2);
                                } else {
                                    break;
                                }
                            }
                        }

                    }
//                    ArrayList<OD_Query> need_prob_queries = new ArrayList<>();
//                    for (OD_Query q : all_queries) {
//                        //check with core points
//                        if (currentTime - q.W < d.arrivalTime
//                                && (currentTime - Constants.W) / q.S * q.S == (currentTime - Constants.W)) {
//                            if (!d.safeInlierQueries.contains(q) && d.total_neighbor(q.R, q.W) < q.k) {
//                                //find core in half r 
//                                CorePoint c = findCorePointInRangeR2(d, q.R, sIdx);
//                                if (c == null || !c.enoughPointsInHalfR(q.R, q.W, q.k)) {
//                                    need_prob_queries.add(q);
//                                }
//                            }
//                        }
//                    }

//                    filterRequestTime += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//                    start = Utils.getCPUTime();
                    //compute k max and r-w map
                    Integer k_max = all_sorted_k.get(all_sorted_k.size() - 1);
                    TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> r_w_map = convertToRWMap(need_prob_queries);

                    HashMap<Double, Integer> r_neighbor_count = new HashMap<>();
                    HashMap<Double, Integer> r_suc_neighbor_count = new HashMap<>();
                    for (Double r : r_w_map.keySet()) {
                        r_neighbor_count.put(r, 0);
                        r_suc_neighbor_count.put(r, 0);
                    }
                    ArrayList<Integer> slide_to_check = new ArrayList<>();
                    for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                        if (r_w_map.keySet().isEmpty()) {
                            break;
                        }
                        Double max_r = r_w_map.lastKey();
                        Double min_r = r_w_map.firstKey();

                        if (d.lastProbRadius.containsKey(idx) && d.lastProbRadius.get(idx) >= max_r
                                && !(idx == last_report_slide && isOverlappingSlide)) {
                            TreeMap<Double, Integer> neighborCountMap = d.get_neighbor_count(d.neighborCount.get(idx),
                                    new ArrayList<>(r_w_map.keySet()));

                            for (Double r : new ArrayList<>(r_w_map.keySet())) {
//                                int n = d.get_neighbor_count(d.neighborCount.get(idx), r);
                                int n = neighborCountMap.get(r);
                                r_neighbor_count.put(r, r_neighbor_count.get(r) + n);
                                if (idx >= d.sIndex) {
                                    r_suc_neighbor_count.put(r, r_suc_neighbor_count.get(r) + n);
                                }

                                checkInlier(r_w_map, r, idx, result, d, r_neighbor_count, r_suc_neighbor_count);

                            }
                        } else {
                            slide_to_check.add(idx);

                        }

                        boolean timeToCheckWindow = false;
                        outer:
                        for (Double r : r_w_map.keySet()) {
                            for (Integer w : r_w_map.get(r).keySet()) {
                                if (currentTime - w == all_slides.get(idx).get(0).arrivalTime - 1) {
                                    timeToCheckWindow = true;
                                    break outer;
                                }
                            }
                        }
                        if (timeToCheckWindow) {
                            while (!slide_to_check.isEmpty() && !r_w_map.isEmpty()) {
                                Integer idx2 = slide_to_check.get(0);
                                probe_slide(d, idx2, max_r, min_r, k_max, r_neighbor_count.get(min_r));
                                TreeMap<Double, Integer> neighborCountMap = d.get_neighbor_count(d.neighborCount.get(idx2),
                                        new ArrayList<>(r_w_map.keySet()));
                                for (Double r : new ArrayList<>(r_w_map.keySet())) {
                                    int n = neighborCountMap.get(r);
                                    r_neighbor_count.put(r, r_neighbor_count.get(r) + n);
                                    if (idx >= d.sIndex) {
                                        r_suc_neighbor_count.put(r, r_suc_neighbor_count.get(r) + n);
                                    }
                                    checkInlier(r_w_map, r, idx, result, d, r_neighbor_count, r_suc_neighbor_count);
                                }
                                slide_to_check.remove(0);
                            }
                            for (Double r : new ArrayList<>( r_w_map.keySet())) {
                                for (Integer w : new ArrayList<>(r_w_map.get(r).keySet())) {
                                    if (currentTime - w == all_slides.get(idx).get(0).arrivalTime - 1) {
                                        for (OD_Query q : (new ArrayList<>(r_w_map.get(r).get(w)))) {
                                            if (r_neighbor_count.get(q.R) < q.k) {
                                                if (!result.containsKey(q)) {
                                                    result.put(q, new ArrayList<>());
                                                }
                                                result.get(q).add(d);
                                            }
                                        }

                                        r_w_map.get(r).remove(w);
                                        if (r_w_map.get(r).isEmpty()) {
                                            r_w_map.remove(r);
                                        }
                                    }
                                }

                            }
                        }

                    }

//                    for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
//                        if (r_w_map.keySet().isEmpty()) {
//                            break;
//                        }
//                        Double max_r = r_w_map.lastKey();
//                        Double min_r = r_w_map.firstKey();
//
//                        if (d.lastProbRadius.containsKey(idx) && d.lastProbRadius.get(idx) >= max_r
//                                && !(idx == last_report_slide && isOverlappingSlide)) {
//
//                        } else {
//                            probe_slide(d, idx, max_r, min_r, k_max, r_neighbor_count.get(min_r));
//                            TreeMap<Double, Integer> neighborCountMap = d.get_neighbor_count(d.neighborCount.get(idx),
//                                    new ArrayList<>(r_w_map.keySet()));
//                            for (Double r : new ArrayList<>(r_w_map.keySet())) {
////                                int n = d.get_neighbor_count(d.neighborCount.get(idx), r);
//                                int n = neighborCountMap.get(r);
//                                r_neighbor_count.put(r, r_neighbor_count.get(r) + n);
//                                if (idx >= d.sIndex) {
//                                    r_suc_neighbor_count.put(r, r_suc_neighbor_count.get(r) + n);
//                                }
//                                checkInlier(r_w_map, r, idx, result, d, r_neighbor_count, r_suc_neighbor_count);
//                            }
//                        }
//
//                    }
//                    findNeighborTime += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
                }
            }

//            System.out.println("Time for filtering requests = " + filterRequestTime);
//            System.out.println("Time for finding neighbors  = " + findNeighborTime);
//            System.out.println("Time for using Core = " + useCoreTime);
//            System.out.println("Time for checkDomiated  = " + checkDominated);
            last_report_slide = newestSlide;
            isOverlappingSlide = false;
            last_report_time = currentTime;
        } else {
            System.out.println("No need reporting!!!");
        }

        return result;
    }

    public TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> convertToRWMap(ArrayList<OD_Query> need_prob_queries) {
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
        return r_w_map;
    }

    public CorePoint findCorePointInRangeR2(C_Data d, Double r, int sIdx) {
        for (CorePoint c : all_core_points_map.get(sIdx).get(r)) {
            if (DistanceFunction.euclideanDistance(d, c) <= r / 2) {
                return c;
            }
        }
        return null;
    }

    public void checkInlier(TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> r_w_map,
            Double r, int idx, HashMap<OD_Query, ArrayList<Data>> result, C_Data d,
            HashMap<Double, Integer> r_neighbor_count, HashMap<Double, Integer> r_suc_neighbor_count) {
        for (Integer w : new ArrayList<>(r_w_map.get(r).keySet())) {
            for (OD_Query q : new ArrayList<>(r_w_map.get(r).get(w))) {
                if (r_neighbor_count.get(r) >= q.k) {
                    r_w_map.get(r).get(w).remove(q);
                    if (r_w_map.get(r).get(w).isEmpty()) {
                        r_w_map.get(r).remove(w);
                        if (r_w_map.get(r).isEmpty()) {
                            r_w_map.remove(r);
                        }
                    }
                    if (r_suc_neighbor_count.get(r) >= q.k) {
                        d.safeInlierQueries.add(q);
                    }
//                    break;
                }

            }
        }
    }

    public boolean isEndOfSlide(int arrivalTime) {

        return arrivalTime >= Constants.W && isStartOfSlide(arrivalTime + 1 - Constants.W);
    }

    public boolean isStartOfSlide(int arrivalTime) {
        if (arrivalTime < 1) {
            return false;
        }
//        for (Integer s : all_s) {
//            if ((arrivalTime - 1) / s * s == (arrivalTime - 1)) {
//                return true;
//            }
//        }
        for (OD_Query q : all_queries) {
            if ((arrivalTime - 1 + q.W - Constants.W >= 0)
                    && (arrivalTime - 1 + q.W - Constants.W) / q.S * q.S == (arrivalTime - 1 + q.W - Constants.W)) {
                return true;
            }
        }
        return false;
    }

    private void probe_slide(C_Data d, int sIdx, Double max_r, Double min_r, int k_max, int count_neigh_min_r) {

        int min_arrival_time = all_slides.get(sIdx).get(0).arrivalTime;
        boolean[] checked = null;
        int case_ = 0;
        ArrayList<Bin> possibleCandidates = new ArrayList<>();
        Double lastR = d.lastProbRadius.get(sIdx);
        if(lastR == null) lastR = -1.0;
        //find close core
        ResultFindCore rf = findCloseCore(d, sIdx, max_r);
        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        if (cores != null) {
            if (distance <= max_r / 2) {
                CorePoint c = cores.get(0);
                //grab close neighbor in range R/2 of c
                possibleCandidates.addAll(c.getDataInRange(0, 1.5 * max_r, sIdx));
//                d.lastProbCore.put(sIdx, c);
            } else if (distance <= max_r) {
                possibleCandidates.addAll(cores.get(0).getDataInRange(0, 2 * max_r, sIdx));
//                d.lastProbCore.put(sIdx, cores.get(0));
            } else if (distance <= max_r * 2) {
                case_ = 1;
                checked = new boolean[all_slides.get(sIdx).size()];
                for (CorePoint c : cores) {
                    possibleCandidates.addAll(c.getDataInRange(0, max_r, sIdx));
                }
//                d.lastProbCore.put(sIdx, null);
            }

        }
//        } 
//        else if (d.lastProbRadius.get(sIdx) < max_r && d.lastProbCore.get(sIdx) != null) {
//
//            CorePoint c = d.lastProbCore.get(sIdx);
//            double distance = DistanceFunction.euclideanDistance(d, c);
//            lastR = d.lastProbRadius.get(sIdx);
//            double lowerbound = lastR / 2;
//            if (distance <= lastR / 2) {
//                lowerbound = lastR / 2;
//            }
//            //else if (distance <= lastR) {
//            //  lowerbound = 2 * lastR;
//            //}
//            if (distance <= max_r / 2) {
//                //grab close neighbor in range R/2 of c
//                possibleCandidates.addAll(c.getDataInRange(lowerbound, 1.5 * max_r, sIdx));
//            } else if (distance <= max_r) {
//                possibleCandidates.addAll(c.getDataInRange(lowerbound, 2 * max_r, sIdx));
//            }
//
//        }

        outer:
        for (Bin ps : possibleCandidates) {
            for (int t = 0; t < ps.data.get(sIdx).size(); t++) {
                C_Data d2 = ps.data.get(sIdx).get(t);
                //for overlapping slide
                if (sIdx == last_report_slide && isOverlappingSlide && d2.arrivalTime <= last_report_time) {
                    continue;
                }

                if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {
                    double dist = DistanceFunction.euclideanDistance(d2, d);
                    if (dist > lastR && dist <= max_r) {
                        Double matchR = getMatchRadius(dist, all_sorted_r);
                        if (!d.neighborCount.containsKey(sIdx)) {
                            d.neighborCount.put(sIdx, new TreeMap<>());
                        }
                        if (!d.neighborCount.get(sIdx).containsKey(matchR)) {
                            d.neighborCount.get(sIdx).put(matchR, 1);
                        } else {
                            d.neighborCount.get(sIdx).put(matchR, d.neighborCount.get(sIdx).get(matchR) + 1);
                        }

                        if (dist <= min_r) {
                            count_neigh_min_r += 1;
                            if (count_neigh_min_r >= k_max) {
                                break outer;
                            }
                        }

                    }

                    if (case_ == 1) {
                        checked[d2.arrivalTime - min_arrival_time] = true;
                    }

                }
            }
        }

        d.lastProbRadius.put(sIdx, max_r);
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

//        if (d.closeCoreMaps_halfR.get(radius) != null
//                && d.closeCoreMaps_halfR.get(radius).supported_slides.contains(slideIndex)) {
//            resultCore = new ArrayList<>();
//            resultCore.add(d.closeCoreMaps_halfR.get(radius));
//            return new ResultFindCore(DistanceFunction.euclideanDistance(d, d.closeCoreMaps_halfR.get(radius)), resultCore);
//        } else if (d.closeCoreMaps_R.get(radius) != null) {
//            CorePoint c = d.closeCoreMaps_R.get(radius);
//
//            if (all_core_points_map.get(slideIndex).get(radius).contains(c)) {
//                resultCore = new ArrayList<>();
//                resultCore.add(c);
//                return new ResultFindCore(DistanceFunction.euclideanDistance(d, c), resultCore);
//            }
//
//        }
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
                    distance_to_cores.add(distance);
                    break;
                } else if (distance <= radius * 2) {
                    inRangeDoubleRCores.add(c);
                    distance_to_cores.add(distance);
                }
            }
        }
        if (!inRangeRCores.isEmpty()) {
            return new ResultFindCore(distance_to_cores.get(0), inRangeRCores);
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
            all_r_w.put(q.R, new TreeMap<>());
        }
        if (!all_r_w.get(q.R).containsKey(q.W)) {
            all_r_w.get(q.R).put(q.W, new ArrayList<>());
        }
        all_r_w.get(q.R).get(q.W).add(q);
//        Collections.sort(all_r_w.get(q.R));
        if (!all_s.contains(q.S)) {
            all_s.add(q.S);
        }
//        for(Integer s: )
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

//        if (all_slides.containsKey(expiredSlideIndex)) {
//            for (C_Data d : all_slides.get(expiredSlideIndex)) {
//                d.reset();
//            }
//        }
        for (int c_idx = all_distinct_core_point.size() - 1; c_idx >= 0; c_idx--) {
            CorePoint c = all_distinct_core_point.get(c_idx);
//            for (Double r : all_sorted_r) {
//                if (r <= c.max_support_r) {
//                    c.discountHalfPoints(r);
//                } else {
//                    break;
//                }
//            }
            for (Bin b : c.all_bins) {
//                if (b.data.containsKey(expiredSlideIndex)) {
                for (Integer sIdx : (new ArrayList<>(b.data.keySet()))) {
                    if (sIdx <= expiredSlideIndex) {
                        b.data.remove(sIdx);
                    }
                }
//                }
            }
            for (Integer sIdx : new ArrayList<>(c.supported_slides)) {
                if (sIdx <= expiredSlideIndex) {
                    c.supported_slides.remove(sIdx);
                }
            }

            if (c.supported_slides.isEmpty()) {
                all_distinct_core_point.remove(c_idx);
            }
        }
        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
//                if (d.neighborCount.containsKey(expiredSlideIndex)) {
                for (Integer slideIndex : new ArrayList<>(d.neighborCount.keySet())) {
                    if (slideIndex <= expiredSlideIndex) {
                        d.neighborCount.remove(slideIndex);
//                        d.lastProbCore.remove(slideIndex);
                        d.lastProbRadius.remove(slideIndex);
                    }
                }

            }
        }
//        if (all_core_points_map.containsKey(expiredSlideIndex)) {
//            all_core_points_map.get(expiredSlideIndex).clear();
//        }
        for (Integer sIdx : new ArrayList<>(all_slides.keySet())) {
            if (sIdx <= expiredSlideIndex) {
                for (C_Data d : all_slides.get(sIdx)) {
                    d.reset();
                }
                all_slides.remove(sIdx);
                all_core_points_map.remove(sIdx);
            }
        }

    }

    private void selectAllCore(int sIdx, ArrayList<Double> all_sorted_r) {
        HashSet<CorePoint> indexed_cores = new HashSet<>();
        ArrayList<C_Data> data = all_slides.get(sIdx);
        ArrayList<CorePoint> current_cores = new ArrayList<>();

        for (int r_idx = all_sorted_r.size() - 1; r_idx >= 0; r_idx--) {
            Double r = all_sorted_r.get(r_idx);

            //for overlapping slide
            if (all_core_points_map.get(sIdx) != null
                    && all_core_points_map.get(sIdx).get(r) != null) {
                for (CorePoint c : all_core_points_map.get(sIdx).get(r)) {
                    if (!current_cores.contains(c)) {
                        current_cores.add(c);
                    }
                }

            }
            for (C_Data d : data) {
                if (d.arrivalTime <= last_report_time) {
                    continue;
                }

                //check with current cores
                boolean foundCore = false;
                for (CorePoint c : current_cores) {
                    double dist = DistanceFunction.euclideanDistance(d, c);
                    if (dist <= r) {
                        //connect d to c
                        foundCore = true;
                        if (!indexed_cores.contains(c)) {
                            c.putDataToBin(d, dist, sIdx);
                        }
                        break;
                    }
                }

                //if d is not linked to a core, find in history
                if (!foundCore) {
                    for (CorePoint c : all_distinct_core_point) {
                        if (c.max_support_r >= r) {
                            double dist = DistanceFunction.euclideanDistance(d, c);
                            if (dist <= r) {
                                foundCore = true;
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
                if (!foundCore) {
                    //creat a core from d
                    CorePoint c = new CorePoint(d);
                    c.creatBins(all_sorted_r.subList(0, r_idx + 1));
                    c.max_support_r = r;
                    current_cores.add(c);
                    c.putDataToBin(d, 0, sIdx);
                    all_distinct_core_point.add(c);
                }
            }

            //index data points around cores
            for (CorePoint c : current_cores) {
                if (!indexed_cores.contains(c)) {
                    boolean[] checked = new boolean[all_slides.get(sIdx).size()];
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
            for (CorePoint c : all_core_points_map.get(sIdx).get(r)) {
                if (!c.supported_slides.contains(sIdx)) {
                    c.supported_slides.add(sIdx);
                }
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
                return 3;
            } else if (o1.R == o2.R) {
                if (o1.k < o2.k) {
                    return 2;
                } else if (o1.k > o2.k) {
                    return -2;
                } else {
                    if (o1.W < o2.W) {
                        return -1;
                    } else if (o1.W > o2.W) {
                        return 1;
                    } else {
                        return 0;
                    }
                }
            } else {
                return -3;
            }
        }

    }

    class CorePoint extends C_Data {

        public ArrayList<Bin> all_bins = new ArrayList<>();
//        public HashMap<Double, Integer> totalHalfRPoints = new HashMap<>();

//        public HashMap<Double, TreeMap<Integer, Integer>> totalHalfRPoints = new HashMap<>();
        public Double max_support_r;
        public HashSet<Integer> supported_slides = new HashSet<>();
//        public HashMap<Double, Integer> lastHalfR = new HashMap<>();

//        public boolean enoughPointsInHalfR(Double r, Integer w, int k) {
//            return computeTotalHalfRPoints(r, w) > k;
//        }
        public boolean enoughPointsInHalfR(Double r, Integer w, int k) {
            int t = 0;
            for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                if (all_slides.get(idx).get(0).arrivalTime >= currentTime + 1 - w) {
                    for (Bin b : all_bins) {
                        if (b.max_val <= r / 2) {
                            if (b.data.containsKey(idx)) {
                                t += b.data.get(idx).size();
                                if (t > k) {
                                    return true;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                } else {
                    return false;
                }

//                if (t > k) {
//                    return true;
//                }
            }

            return false;
        }

        public int computeTotalHalfRPoints(Double r, int w) {
            int t = 0;
            for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                if (all_slides.get(idx).get(0).arrivalTime >= currentTime + 1 - w) {
                    for (Bin b : all_bins) {
                        if (b.max_val <= r / 2) {
                            if (b.data.containsKey(idx)) {
                                t += b.data.get(idx).size();
                            }
                        } else {
                            break;
                        }
                    }
                }

                if (t > all_sorted_k.get(all_sorted_k.size() - 1)) {
                    break;
                }
            }

            return t;
        }

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
        //public HashMap<Double, CorePoint> closeCoreMaps_halfR = new HashMap<>();
        //public HashMap<Double, CorePoint> closeCoreMaps_R = new HashMap<>();
        public HashMap<Integer, Double> lastProbRadius = new HashMap<>();
//        public HashMap<Integer, CorePoint> lastProbCore = new HashMap<>();

        public int sIndex = -1;
        public HashSet<OD_Query> safeInlierQueries = new HashSet<>();

//        public void connectCore(CorePoint c, double distance, List<Double> all_r) {
//            for (int r_idx = all_r.size() - 1; r_idx >= 0; r_idx--) {
//                Double r = all_r.get(r_idx);
//                if (distance <= r) {
//                    if (distance <= r / 2) {
//                        closeCoreMaps_halfR.put(r, c);
//                    } else {
//                        closeCoreMaps_R.put(r, c);
//                    }
//
//                } else {
//                    break;
//                }
//            }
//        }
        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.sIndex = (arrivalTime - 1) / Constants.slide;
        }

        public C_Data() {

        }

        public int total_neighbor(Double r, int w) {
            int n = 0;
            for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                if (lastProbRadius.containsKey(idx) && lastProbRadius.get(idx) >= r
                        && !(idx == last_report_slide && isOverlappingSlide)
                        && all_slides.get(idx).get(0).arrivalTime >= currentTime - w) {
                    n += get_neighbor_count(neighborCount.get(idx), r);
                }

            }
            return n;
        }

        public TreeMap<Double, Integer> total_suc_neighbor_map(TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> r_w_map) {
            TreeMap<Double, Integer> countMap = new TreeMap<>();
            for (Double r : r_w_map.keySet()) {
                countMap.put(r, 0);
            }
            for (int idx = newestSlide; idx >= sIndex; idx--) {
                TreeMap<Double, Integer> r_count_map = get_neighbor_count(neighborCount.get(idx), new ArrayList<>(r_w_map.keySet()));
                for (Double r : r_count_map.keySet()) {
                    countMap.put(r, countMap.get(r) + r_count_map.get(r));
                }
            }
            return countMap;
        }

        public TreeMap<Double, TreeMap<Integer, Integer>> total_neighborMap(TreeMap<Double, TreeMap<Integer, ArrayList<OD_Query>>> r_w_map) {
            TreeMap<Double, TreeMap<Integer, Integer>> countMap = new TreeMap<>();
            for (Double r : r_w_map.keySet()) {
                TreeMap<Integer, Integer> w_map = new TreeMap<>();
                for (Integer w : r_w_map.get(r).keySet()) {
                    w_map.put(w, 0);
                }
                countMap.put(r, w_map);
            }
            if (!r_w_map.isEmpty()) {
                for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                    TreeMap<Double, Integer> r_count_map = get_neighbor_count(neighborCount.get(idx), new ArrayList<>(r_w_map.keySet()));
                    for (Double r : r_count_map.keySet()) {
                        for (Integer w : countMap.get(r).descendingKeySet()) {
                            if (all_slides.get(idx).get(0).arrivalTime >= currentTime - w) {
                                countMap.get(r).put(w, countMap.get(r).get(w) + r_count_map.get(r));
                            } else {
                                break;
                            }

                        }
                    }
                }
            }
            return countMap;
        }

        public boolean has_enough_neighbor(Double r, int w, int k) {
            int n = 0;
            for (int idx = newestSlide; idx > expiredSlideIndex; idx--) {
                if (lastProbRadius.containsKey(idx) && lastProbRadius.get(idx) >= r
                        && !(idx == last_report_slide && isOverlappingSlide)
                        && all_slides.get(idx).get(0).arrivalTime >= currentTime - w) {
                    n += get_neighbor_count(neighborCount.get(idx), r);

                }
                if (n >= k) {
                    return true;
                }
                if (all_slides.get(idx).get(0).arrivalTime < currentTime - w) {
                    return false;
                }

            }
            return false;
        }

        private Integer get_neighbor_count(TreeMap<Double, Integer> neighbor_count_map, Double r) {
            if (neighbor_count_map == null) {
                return 0;
            }
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

        private TreeMap<Double, Integer> get_neighbor_count(TreeMap<Double, Integer> neighbor_count_map, ArrayList<Double> r_list) {
            TreeMap<Double, Integer> result = new TreeMap<>();
            for (Double r : r_list) {
                result.put(r, 0);
            }
            if (neighbor_count_map == null) {
                return result;
            }
            int count = 0;
            int cur_idx = 0;
            for (Double ri : neighbor_count_map.keySet()) {
                count += neighbor_count_map.get(ri);
                if (cur_idx < r_list.size()) {
                    if (Objects.equals(ri, r_list.get(cur_idx))) {
                        result.put(r_list.get(cur_idx), count);
                        cur_idx += 1;
                    }
                } else {
                    break;
                }

            }
            //update the rest of r_list
            for (int i = cur_idx; i < r_list.size(); i++) {
                result.put(r_list.get(i), count);
            }
            return result;
        }

        public HashMap<Double, Integer> neighborCounts(ArrayList<Double> r_list) {
            HashMap<Double, Integer> result = new HashMap<>();
            for (Integer sIdx : neighborCount.keySet()) {
                TreeMap<Double, Integer> cur_neigh = get_neighbor_count(neighborCount.get(sIdx), r_list);
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

        private void reset() {
//            this.lastProbCore.clear();
            this.lastProbRadius.clear();
            this.neighborCount.clear();
            this.safeInlierQueries.clear();
        }

    }
}
