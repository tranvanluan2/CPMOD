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
import mtree.tests.Data;
import mtree.utils.Constants;

/**
 *
 * @author luan
 */
public class MCSKy {

    public static ArrayList<Micro_Cluster> microclusters = new ArrayList<>();
    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();
    public ArrayList<Double> all_r;
    public ArrayList<Integer> all_k;
    public HashMap<Integer, ArrayList<OD_Query>> k_query_map = new HashMap<>();
    public static int minW = Integer.MAX_VALUE;

    public void add_query(OD_Query q) {
        all_queries.add(q);

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
        //update min window 
        if (minW > q.W) {
            minW = q.W;
        }

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

    public void create_hashmap_k_queries() {
        for (OD_Query q : all_queries) {
            if (k_query_map.containsKey(q.k)) {

                k_query_map.get(q.k).add(q);
            } else {
                ArrayList<OD_Query> queries = new ArrayList<>();
                queries.add(q);
                k_query_map.put(q.k, queries);
            }
        }
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

    public static HashMap<Integer, ArrayList<C_Data>> micro_clusters = new HashMap<>();

    public HashMap<OD_Query, ArrayList<Data>> slide_process(ArrayList<Data> data, int _currentTime) {

        HashMap<OD_Query, ArrayList<Data>> result = new HashMap<>();
        for (OD_Query q : all_queries) {
            result.put(q, new ArrayList<>());
        }
        currentTime = _currentTime;
        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
        int[] slide_to_process;
        if (currentTime == Constants.W) {
            slide_to_process = new int[Constants.W / Constants.slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
            get_unique_r();
            create_hashmap_k_queries();

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

        expiredSlideIndex = (currentTime - 1) / Constants.slide - Constants.W / Constants.slide;
        processExpiredData(expiredSlideIndex);
        int newestSlide = (currentTime - 1) / Constants.slide;
        for (Integer sidx : slide_to_process) {
            for (C_Data d : all_slides.get(sidx)) {
                if (d.sIndex > newestSlide - minW / Constants.slide) //find micro-cluster to add 
                {
                    for (Micro_Cluster c : microclusters) {
                        double distance = DistanceFunction.euclideanDistance(d, c.center);
                        if (distance <= all_r.get(0) / 2.0) {
                            c.members.add(d);
                            d.mc = c;
                            break;
                        }
                    }
                }
            }
        }
        //
        for (int sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                if (d.mc == null) {
                    ArrayList<OD_Query> need_evaluate_queries = new ArrayList<>();
                    for (OD_Query q : all_queries) {
                        if (!d.safe_inlier_queries.contains(q) && d.arrivalTime > currentTime - q.W) {
                            need_evaluate_queries.add(q);
                        }
                    }
                    if (need_evaluate_queries.size() > 0) {
                        HashSet<OD_Query> inlier_queries = create_skyband(d, need_evaluate_queries);
                        for (OD_Query q : need_evaluate_queries) {
                            if ((currentTime - Constants.W) / q.S *q.S ==(currentTime - Constants.W) 
                                    && !inlier_queries.contains(q)) {
                                result.get(q).add(d);
                            }
                        }
                    }
                }
            }

        }
        return result;

    }

    private void processExpiredData(int expiredSlideIndex) {
        all_slides.remove(expiredSlideIndex);

        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                if (d.mc == null && d.sky != null) {
                    for (Bucket b : d.sky.buckets) {
                        ArrayList<C_Data> need_remove = new ArrayList<>();
                        for (int i = 0; i < b.data.size(); i++) {
                            C_Data d2 = b.data.get(i);
                            if (d2.sIndex <= expiredSlideIndex) {
                                need_remove.add(d2);
                            }
                        }
                        for (C_Data d2 : need_remove) {
                            b.data.remove(d2);
                        }

                    }
                }

            }
        }

        int expiredSlideForMC = (currentTime - 1) / Constants.slide - minW/Constants.slide;
        ArrayList<Micro_Cluster> to_remove = new ArrayList<>();
        for (Micro_Cluster mc : microclusters) {

            //remove expired members
            ArrayList<C_Data> member_remove = new ArrayList<>();
            for (C_Data d : mc.members) {
                if (d.sIndex <= expiredSlideForMC) {
                    member_remove.add(d);
                }
            }
            for (C_Data d : member_remove) {
                mc.members.remove(d);
            }

            if (mc.members.size() <= all_k.get(0)) {
                //disperse cluster
                to_remove.add(mc);
                for (C_Data d : mc.members) {
                    d.mc = null;
                    d.reset();
                }
            }
        }
        for (Micro_Cluster mc : to_remove) {
            microclusters.remove(mc);
        }
    }

    private int normalize(double distance) {
        int dist = -1;
        for (int i = 0; i < all_r.size(); i++) {
            if (distance <= all_r.get(i)) {
                dist = i;
                break;
            }
        }
        return dist;
    }

    class Micro_Cluster {

        public C_Data center;
        public ArrayList<C_Data> members = new ArrayList<>();
    }

    class C_Data extends Data {

        public int sIndex = -1;
        public Skyband sky;
        public HashSet<OD_Query> safe_inlier_queries = new HashSet<>();
        public Micro_Cluster mc;

//        public HashSet<C_Data> succeeding_neighbors = new HashSet<>();
//
//        public HashSet<C_Data> preceding_neighbors = new HashSet<>();
        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;
        }

        public ArrayList<C_Data> getSortedSkybandByTime() {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Bucket b : sky.buckets) {
                result.addAll(b.data);
            }
            Collections.sort(result, new Comparator<C_Data>() {
                @Override
                public int compare(C_Data o1, C_Data o2) {
                    if (o1.arrivalTime > o2.arrivalTime) {
                        return -1;
                    } else if (o1.arrivalTime == o2.arrivalTime) {
                        return 0;
                    } else {
                        return 1;
                    }
                }
            });
            return result;
        }

        private boolean skyEvaluate(C_Data d, double distance,
                HashSet<OD_Query> inlier_queries, ArrayList<OD_Query> need_evaluate_queries) {
            int normalized_distance = normalize(distance);
            int count = 0;
            Bucket bm = null;
            for (Bucket b : sky.buckets) {
                if (b.normalized_distance <= normalized_distance) {
                    count += b.data.size();
                    if (b.normalized_distance == normalized_distance) {
                        bm = b;
                    }
                }

            }
            if (bm != null) {

                double maxrn = need_evaluate_queries.get(need_evaluate_queries.size() - 1).R;

                int max_k = -1;
                for (OD_Query q : need_evaluate_queries) {
                    if (max_k < q.k) {
                        max_k = q.k;
                    }
                }
                if (count < max_k && distance <= maxrn) {
                    bm.data.add(d);

                    //check with k of queries to decide inlier
                    if (k_query_map.containsKey(count + 1)) {
                        for (OD_Query q : k_query_map.get(count + 1)) {
                            if (d.arrivalTime > currentTime - q.W && 
                                    q.R >= distance && need_evaluate_queries.contains(q)) {
                                inlier_queries.add(q);
                                if (d.sIndex >= this.sIndex) {
                                    safe_inlier_queries.add(q);
                                }
                            }

                        }
                    }
                    return true;
                }
                else return false;
            }
            return true;
        }

        private void reset() {

            //create buckets or reset
            if (sky == null) {
                sky = new Skyband();
                sky.buckets = new ArrayList<>();
                for (int i = 0; i < all_r.size(); i++) {
                    Bucket b = new Bucket(i);
                    sky.buckets.add(b);
                }
            } else {
                for (Bucket b : sky.buckets) {
                    b.data.clear();
                }
            }
        }

    }

    private HashSet<OD_Query> create_skyband(C_Data d, ArrayList<OD_Query> need_evaluate_queries) {
        ArrayList<C_Data> input = new ArrayList<>();
        if (d.sky == null) {
            //input is the current window
            expiredSlideIndex = (currentTime - 1) / Constants.slide - Constants.W / Constants.slide;
            int newestSlide = (currentTime - 1) / Constants.slide;
            for (int sidx = newestSlide; sidx > expiredSlideIndex; sidx--) {
                input.addAll(all_slides.get(sidx));
            }
        } else {
            //input is the new slide + existing skyband  
            int newestSlide = (currentTime - 1) / Constants.slide;
            input.addAll(all_slides.get(newestSlide));
            input.addAll(d.getSortedSkybandByTime());
        }
        d.reset();
        HashSet<OD_Query> inlier_queries = new HashSet<>();
        int newestSlide = (currentTime - 1) / Constants.slide;
        ArrayList<C_Data> to_form_cluster = new ArrayList<>();
        for (C_Data d2 : input) {
            double distance = DistanceFunction.euclideanDistance(d, d2);
            //check to form cluster
            if (d.sIndex > newestSlide - minW/Constants.slide
                    && d2.sIndex > newestSlide - minW/Constants.slide
                    && d2.mc == null && distance <= all_r.get(0) / 2) {
                to_form_cluster.add(d2);
                if (to_form_cluster.size() > all_k.get(all_k.size() - 1)) {
                    //form new cluster
                    Micro_Cluster mc = new Micro_Cluster();
                    mc.center = d;
                    mc.members.add(d);
                    mc.members.addAll(to_form_cluster);
                    microclusters.add(mc);
                    for (C_Data m : mc.members) {
                        m.mc = mc;
                        m.reset();
                    }
                    inlier_queries.addAll(need_evaluate_queries);
                    break;
                }
            }
            if (d.skyEvaluate(d2, distance, inlier_queries, need_evaluate_queries) == true) {

            } else {
                if (distance <= need_evaluate_queries.get(0).R) {
                    break;
                }
            }
        }
        return inlier_queries;
    }

    class Skyband {

        public ArrayList<Bucket> buckets;

    }

    class Bucket {

        int normalized_distance = -1;
        ArrayList<C_Data> data = new ArrayList<>();

        public Bucket(int normalized_dist) {
            normalized_distance = normalized_dist;
        }
    }

}
