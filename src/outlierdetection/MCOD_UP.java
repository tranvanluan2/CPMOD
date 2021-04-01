/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import be.tarsos.lsh.Index;
import be.tarsos.lsh.LSH;
import be.tarsos.lsh.Vector;
import be.tarsos.lsh.families.EuclidianHashFamily;
import be.tarsos.lsh.families.HashFamily;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class MCOD_UP {

    class ULObject extends Data {

        public boolean isOutlier = false;
        public boolean isSafe;
        public int numSucEvidence = 0;
        public int numPreEvidence = 0;

        public double mean;
        public double std;
        public double[] a;
        public int last_probed_slide_left = -1;
        public int last_probed_slide_right = -1;
        public int last_up_probed_slide_left = -1;
        public int last_up_probed_slide_right = -1;

        public MicroCluster cluster;
        public Slide slide;
        public HashMap<Integer, Integer> neighborCount = new HashMap<>();

        public ULObject(Data d, int currentTime) {

            super();

            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.mean = 0;
            for (int i = 0; i < this.values.length; i++) {
                this.mean += this.values[i];
            }
            this.mean = this.mean / this.values.length;

            this.std = 0;
            for (int i = 0; i < this.values.length; i++) {
                this.std += (this.values[i] - this.mean) * (this.values[i] - this.mean);
            }
            this.std = this.std / this.values.length;
            this.std = Math.sqrt(this.std);

            this.a = new double[this.values.length];
            for (int i = 0; i < this.values.length; i++) {
                this.a[i] = this.values[i] - this.mean;
            }
        }
    }

    class Slide {

        public List<ULObject> points = new ArrayList<>();
        public int id;

        public Window window;

        public LSH lsh;
        HashFamily family;
        public HashMap< ULObject, Vector> ulv = new HashMap<>();

        public HashMap<Vector, ULObject> v_ul = new HashMap<>();

        public Slide(ArrayList<Data> data, int start_time) {
            for (int i = 0; i < data.size(); i++) {
                ULObject obj = new ULObject(data.get(i), start_time + i);
                obj.slide = this;
                points.add(obj);

            }

            ArrayList<Vector> datas = new ArrayList<>();
            points.forEach((d) -> {
                Vector v = new Vector(d);
                datas.add(v);
                ulv.put(d, v);
                v_ul.put(v, d);
            });

            if ((int) (1 * Constants.R) == 0) {
                this.family = new EuclidianHashFamily(4, Constants.dimensions);
            } else {
                this.family = new EuclidianHashFamily((int) (10 * Constants.R), Constants.dimensions);
            }

//            this.family = new EuclidianHashFamily(10000, Constants.dimensions);
            if (lsh == null) {
                lsh = new LSH(datas, family);
                lsh.index = new Index(family, 50, 100);

            }

//            datas.forEach((d) -> {
//                lsh.index.index(d);
//            });
        }

    }

    class Window {

        public ArrayList<Slide> slides = new ArrayList<>();

        public Window() {
        }

        public Slide getNewestSlide() {
            if (slides.isEmpty()) {
                return null;
            }
            return slides.get(slides.size() - 1);
        }

        public Slide getExpiringSlide() {
            if (cur_slide_id <= Constants.W / Constants.slide) {
                return null;
            } else {
                return slides.get(0);
            }
        }

        public void addNewSlide(Slide s) {
            Slide newestSlide = getNewestSlide();
            if (newestSlide != null) {
                s.id = newestSlide.id + 1;
            } else {
                s.id = 0;
            }
            s.window = this;
            slides.add(s);
        }

    }

    class MicroCluster {

        public ArrayList<ULObject> members;
        public ULObject center;

        public MicroCluster(ULObject center, ArrayList<ULObject> members) {
            this.center = center;
            this.members = members;
            this.members.add(this.center);
            members.forEach((o) -> {
                o.cluster = this;
            });
        }

        public boolean isDispersed() {
            return this.members.size() < Constants.k + 1;
        }

    }

    //only use for EU distance
    int level_distance = 40;

    public double UB_distance(ULObject x, ULObject y) {
        double distance = 0;
        int d = x.dimensions();

        distance += d * ((x.mean - y.mean) * (x.mean - y.mean) + (x.std + y.std) * (x.std + y.std));
        double epsilon = 0.000000001;

//        System.out.println("x.std = "+x.std+", ");
        for (int i = d - 1; i >= d - level_distance; i--) {
            distance -= x.std * y.std * (y.a[i] / (y.std + epsilon) + x.a[i] / (x.std + epsilon)) * (y.a[i] / (y.std + epsilon) + x.a[i] / (x.std + epsilon));
//            System.out.println("i = " + i +", distance = "+Math.sqrt(distance));

        }
        return Math.sqrt(distance);

    }

    public double LB_distance(ULObject x, ULObject y) {
        double distance = 0;
        int d = x.dimensions();

        distance += d * ((x.mean - y.mean) * (x.mean - y.mean) + (x.std - y.std) * (x.std - y.std));

        double epsilon = 0.000000001;
        for (int i = 0; i < level_distance; i++) {
            distance += x.std * y.std * (y.a[i] / (y.std + epsilon) - x.a[i] / (x.std + epsilon)) * (y.a[i] / (y.std + epsilon) - x.a[i] / (x.std + epsilon));
        }
        return Math.sqrt(distance);
    }

    public double EU_ditance(ULObject x, ULObject y) {
        double distance = 0;
        int d = x.dimensions();
        for (int i = 0; i < d; i++) {
            distance += (x.values[i] - y.values[i]) * (x.values[i] - y.values[i]);
        }
        return Math.sqrt(distance);

    }

    public static ArrayList<MicroCluster> microClusters = new ArrayList<>();
    public static HashSet<ULObject> PDList = new HashSet<>();
    public Window window = new Window();
    public static int cur_slide_id = 0;
    public static int start_time = 0;
    public static HashMap<Integer, ArrayList<ULObject>> triggerList = new HashMap<>();
    public static int use_ub = 0;
    public static int added_to_cluster = 0;
    public static HashSet<ULObject> checked_by_clusters = new HashSet<>();

    public void findNeighborNaive() {
        window.slides.forEach((s) -> {
            s.points.stream().map((o) -> {
                window.slides.forEach((s2) -> {
                    s2.points.stream().filter((o2) -> (o2.arrivalTime != o.arrivalTime && EU_ditance(o, o2) <= Constants.R)).forEachOrdered((_item) -> {
                        if (s2.id >= s.id) {
                            o.numSucEvidence += 1;
                        } else {
                            o.numPreEvidence += 1;

                        }
                    });
                });
                return o;
            }).forEachOrdered((o) -> {
                PDList.add(o);
            });
        });
    }

    public void findNeighborFirstWindow() {
        for (Slide s : window.slides) {
            for (ULObject o : s.points) {

                if (o.cluster != null) {
                    continue;
                }

                //find cluster for o
                for (MicroCluster c : microClusters) {
                    if (o.cluster == null && EU_ditance(o, c.center) <= Constants.R / 2) {
                        c.members.add(o);
                        o.cluster = c;
                        break;
                    }
                }
                if (o.cluster != null) {
                    continue;
                }
                //form cluster for o
                ArrayList<ULObject> closeNeighbors = new ArrayList<>();
                for (Slide s2 : window.slides) {
                    for (ULObject o2 : s2.points) {
                        if (o2.arrivalTime != o.arrivalTime) {
                            if (o2.cluster == null && EU_ditance(o, o2) <= Constants.R / 2) {
                                closeNeighbors.add(o2);
                                if (closeNeighbors.size() >= Constants.k) {
                                    MicroCluster c = new MicroCluster(o, closeNeighbors);
                                    microClusters.add(c);
                                    break;
                                }
                            }
                        }

                    }
                }

                if (o.cluster != null) {
                    continue;
                }
                //find neighbor for o and add o to PD List
                for (int i = window.slides.size() - 1; i >= 0; i--) {
                    Slide s2 = window.slides.get(i);
                    for (ULObject o2 : s2.points) {
                        //test
//                        System.out.println("UB distance = " + UB_distance(o, o2));
//                        System.out.println("LB distance = " + LB_distance(o, o2));
//                        System.out.println("Exact distance = " + EU_ditance(o, o2));
                        //end test
                        if (o2.arrivalTime != o.arrivalTime && EU_ditance(o, o2) <= Constants.R) {
//                    s2.points.stream().filter((o2) -> (o2.arrivalTime != o.arrivalTime && EU_ditance(o, o2) <= Constants.R)).forEachOrdered((_item) -> {
                            if (s2.id >= s.id) {
                                o.numSucEvidence += 1;
                            } else {
                                o.numPreEvidence += 1;
                                if (o.neighborCount.containsKey(s2.id)) {
                                    o.neighborCount.put(s2.id, o.neighborCount.get(s2.id) + 1);
                                } else {
                                    o.neighborCount.put(s2.id, 1);
                                }
                                if (triggerList.containsKey(s2.id)) {
                                    ArrayList<ULObject> points = triggerList.get(s2.id);
                                    points.add(o);
                                } else {
                                    ArrayList<ULObject> points = new ArrayList<>();
                                    points.add(o);
                                    triggerList.put(s2.id, points);
                                }
                            }
                        }
//                    });
                        if (o.numPreEvidence + o.numSucEvidence >= Constants.k) {
                            break;
                        }
                    }

                }
                if (o.cluster == null) {
                    PDList.add(o);
                }
            }

        }

        //Remove data points that are moved to cluster from PD list
        ArrayList<ULObject> remove_list = new ArrayList<>();
        PDList.stream().filter((o) -> (o.cluster != null)).forEachOrdered((o) -> {
            remove_list.add(o);
        });
        remove_list.forEach((o) -> {
            PDList.remove(o);
        });

    }

    Comparator<ULObject> ulobjectComparator = (ULObject s1, ULObject s2) -> s1.arrivalTime - s2.arrivalTime;

    private void findNeighbor(ULObject o) {

        PriorityQueue<ULObject> neighbors = new PriorityQueue<>(ulobjectComparator);
        ArrayList<ULObject> closeNeighbors = new ArrayList<>();
        HashSet<ULObject> ub_checked = new HashSet<>();

        int numSlide = window.slides.size();

        //Find neighbors using LSH + upper bound 
//        for (int i = numSlide - 1; i >= 0; i--) {
//            int count_success = 0;
//            Slide s = window.slides.get(i);
//            List<Vector> vectors = s.lsh.query(o.slide.ulv.get(o), Constants.k);
////            System.out.println("Size of results = " + vectors.size());
//            for (Vector v : vectors) {
//                ULObject o2 = s.v_ul.get(v);
//                if (o2.arrivalTime != o.arrivalTime) {
//                    //compute upper distance
//                    double ub_d = UB_distance(o, o2);
////                    double exact_distance = EU_ditance(o, o2);
//                    if (ub_d <= Constants.R) {
////                        System.out.println("Good!");
//                        count_success += 1;
//                        neighbors.add(o2);
//                        if (!o2.isSafe && o2.slide.id < o.slide.id) {
//                            o2.numSucEvidence += 1;
//                            if (o2.numSucEvidence >= Constants.k) {
//                                o2.isSafe = true;
//                            }
//                        }
//                        ub_checked.add(o2);
//                        if (ub_d <= Constants.R / 2) {
//                            closeNeighbors.add(o2);
//                        }
//                    } 
////                    else {
//////                        System.out.println("ubd distance = " + ub_d);
////                        break;
////                    }
//                }
//                if (closeNeighbors.size() >= Constants.k) {
//                    MicroCluster c = new MicroCluster(o, closeNeighbors);
//                    microClusters.add(c);
////                    System.out.println("AAAAAAAAA");
//                    break;
//                }
//            }
//            if (vectors.size() > 0) {
////                System.out.println("Count success = " + (count_success * 1.0 / vectors.size()));
//            }
//
//        }
        //Find neighbors for o by going through slides and then add o to PD List
        //find neighbors using upper bound distance 
        if (neighbors.size() < Constants.k) {
            for (int i = numSlide - 1; i >= 0; i--) {
                Slide s = window.slides.get(i);

                for (ULObject o2 : s.points) {
                    if (o2.arrivalTime != o.arrivalTime && !checked_by_clusters.contains(o2)) {

                        double ub_d = UB_distance(o, o2);
                        if (ub_d <= Constants.R) {
                            ub_checked.add(o2);
                            neighbors.add(o2);
                            if (!o2.isSafe && o2.slide.id < o.slide.id) {
                                o2.numSucEvidence += 1;
                                if (o2.numSucEvidence >= Constants.k) {
                                    o2.isSafe = true;
                                }
                            }

//                        if (ub_d <= Constants.R / 2) {
//                            closeNeighbors.add(o2);
//                        }
                        } else {

                        }
                    }

//                if (closeNeighbors.size() >= Constants.k) {
//                    MicroCluster c = new MicroCluster(o, closeNeighbors);
//                    microClusters.add(c);
////                    System.out.println("AAAAAAAAA");
//                    break;
//                }
//                if (neighbors.size() >= 1.5 * Constants.k) {
//                    break;
//
//                }
//                if (neighbors.size() >= Constants.k) {
//                    //stop probing
//                    use_ub += 1;
//
//                    break;
//                }
                }
//            o.last_up_probed_slide_left = s.id;
                if (neighbors.size() >= Constants.k) {
                    break;

                }

            }
        }

        if (neighbors.size() < Constants.k) //find neighbors using lb+exact distance
        {
            for (int i = numSlide - 1; i >= 0; i--) {
                Slide s = window.slides.get(i);
                int totalPoint = 0;
                int filterByLB = 0;
                for (ULObject o2 : s.points) {

                    if (o2.arrivalTime != o.arrivalTime && !ub_checked.contains(o2)
                            && !checked_by_clusters.contains(o2)) {
                        totalPoint += 1;
                        //compute LB  distance
                        double lb_d = LB_distance(o, o2);
                        if (lb_d <= Constants.R) {
                            double exact_d = EU_ditance(o, o2);
                            if (exact_d <= Constants.R) {
                                //o2 is a neighbor of o
                                neighbors.add(o2);
                                if (!o2.isSafe && o2.slide.id < o.slide.id) {
                                    o2.numSucEvidence += 1;
                                    if (o2.numSucEvidence >= Constants.k) {
                                        o2.isSafe = true;
                                    }
                                }

                                if (exact_d <= Constants.R / 2) {
                                    closeNeighbors.add(o2);
                                }
                            }
                        } else {
                            filterByLB += 1;
                        }

                    }

                    if (closeNeighbors.size() >= Constants.k) {
                        MicroCluster c = new MicroCluster(o, closeNeighbors);
                        microClusters.add(c);
                        break;
                    }

//                    if (neighbors.size() >= 1.5 * Constants.k) {
//                        break;
//                    }
//                    if (neighbors.size() >= Constants.k) {
//                        break;
//
//                    }
                }
//                System.out.println("Filterd by LB = "+ filterByLB*1.0/totalPoint);

                o.last_probed_slide_left = s.id;

                //check if found enough  neighbors
                if (neighbors.size() >= Constants.k) {
                    break;

                }
            }
        }

        if (o.cluster == null) {
            //update neighbor list for o
            o.neighborCount = new HashMap<>();
            while (neighbors.size() > Constants.k) {
                neighbors.poll();
            }
            neighbors.stream().map((o2) -> o2.slide.id).forEachOrdered((slide_id) -> {
                if (slide_id >= o.slide.id) {
                    o.numSucEvidence += 1;
                } else {
                    o.numPreEvidence += 1;
                    if (o.neighborCount.containsKey(slide_id)) {
                        o.neighborCount.put(slide_id, o.neighborCount.get(slide_id) + 1);
                    } else {
                        o.neighborCount.put(slide_id, 1);
                    }
                    if (triggerList.containsKey(slide_id)) {
                        ArrayList<ULObject> points = triggerList.get(slide_id);
                        points.add(o);
                    } else {
                        ArrayList<ULObject> points = new ArrayList<>();
                        points.add(o);
                        triggerList.put(slide_id, points);
                    }
                }
            });
        }

        if (o.numSucEvidence >= Constants.k) {
            o.isSafe = true;
        }

        if (o.cluster == null) {
            PDList.add(o);
        }

    }

    private void reprob(ULObject o) {

        //reset neighbor info
        o.numPreEvidence = 0;
        o.numSucEvidence = 0;
        o.neighborCount.clear();

        HashSet<ULObject> ub_checked = new HashSet<>();
        int num_slide = window.slides.size();

        //probing current and succeedinng slide
        for (int i = num_slide - 1; i >= 0; i--) {
            Slide s = window.slides.get(i);
            if (s.id >= o.slide.id) {
                //using upper distance  first
                s.points.stream().filter((o2) -> (o2.arrivalTime != o.arrivalTime)).forEachOrdered((o2) -> {
                    double ub_d = UB_distance(o, o2);
                    if (ub_d <= Constants.R) {
                        ub_checked.add(o2);
                        o.numSucEvidence += 1;
                    }
                });
            } else {
                s.points.forEach((o2) -> {
                    double ub_d = UB_distance(o, o2);
                    if (ub_d <= Constants.R) {
                        ub_checked.add(o2);
                        o.numPreEvidence += 1;
                        if (o.neighborCount.containsKey(s.id)) {
                            o.neighborCount.put(s.id, o.neighborCount.get(s.id) + 1);

                        } else {
                            o.neighborCount.put(s.id, 1);
                        }

                        if (triggerList.containsKey(s.id)) {
                            ArrayList<ULObject> points = triggerList.get(s.id);
                            points.add(o);
                        } else {
                            ArrayList<ULObject> points = new ArrayList<>();
                            points.add(o);
                            triggerList.put(s.id, points);
                        }
                    }
                });
            }

            if (o.numSucEvidence >= Constants.k) {
                o.isSafe = true;
            }

            if (o.numPreEvidence + o.numSucEvidence >= Constants.k) {
                break;
            }

        }

        if (o.numPreEvidence + o.numSucEvidence < Constants.k) {
            //using lb + exact 
            for (int i = num_slide - 1; i >= 0; i--) {
                Slide s = window.slides.get(i);
                if (s.id >= o.slide.id) {
                    s.points.stream().filter((o2) -> (o2.arrivalTime != o.arrivalTime && !ub_checked.contains(o2))).forEachOrdered((o2) -> {
                        double lb_d = LB_distance(o, o2);
                        if (lb_d <= Constants.R) {
                            double eu_d = EU_ditance(o, o2);
                            if (eu_d <= Constants.R) {
                                o.numSucEvidence += 1;
                            }
                        }
                    });
                } else {
                    s.points.stream().filter((o2) -> (o2.arrivalTime != o.arrivalTime && !ub_checked.contains(o2))).forEachOrdered((o2) -> {
                        double lb_d = LB_distance(o, o2);
                        if (lb_d <= Constants.R) {
                            double eu_d = EU_ditance(o, o2);
                            if (eu_d <= Constants.R) {
                                o.numPreEvidence += 1;
                                if (o.neighborCount.containsKey(s.id)) {
                                    o.neighborCount.put(s.id, o.neighborCount.get(s.id) + 1);
                                } else {
                                    o.neighborCount.put(s.id, 1);
                                }

                                if (triggerList.containsKey(s.id)) {
                                    ArrayList<ULObject> points = triggerList.get(s.id);
                                    points.add(o);
                                } else {
                                    ArrayList<ULObject> points = new ArrayList<>();
                                    points.add(o);
                                    triggerList.put(s.id, points);
                                }
                            }
                        }
                    });
                }
                if (o.numPreEvidence + o.numSucEvidence >= Constants.k) {
                    break;
                }
            }
        }

    }

    public ArrayList<Data> detectOutlier(ArrayList<Data> data) {

        long start = Utils.getCPUTime();
        int count = 0;
        ArrayList<Data> outlierList = new ArrayList<>();
        if (data.size() == Constants.W) {
            // split into slides
            int numSlide = (int) Math.ceil(Constants.W * 1.0 / Constants.slide);

            for (int i = 0; i < numSlide; i++) {

                ArrayList<Data> d = new ArrayList<>();
                for (int j = 0; j < Constants.slide; j++) {
                    if (i * Constants.slide + j < data.size()) {
                        d.add(data.get(i * Constants.slide + j));
                    }
                }
                Slide s = new Slide(d, start_time);
                start_time += Constants.slide;
                s.id = cur_slide_id;
                cur_slide_id += 1;
                window.addNewSlide(s);

            }

        } else if (data.size() == Constants.slide) {
            // add this slide to window
            Slide s = new Slide(data, start_time);
            start_time += Constants.slide;
            s.id = cur_slide_id;
            cur_slide_id += 1;
            window.addNewSlide(s);
        }

//        System.out.println("Adding new slide time = " + (Utils.getCPUTime() - start));
        start = Utils.getCPUTime();
        //Remove expire slide
        Slide expiredSlide = window.getExpiringSlide();
        if (expiredSlide != null) {

            //Remove expired data from micro cluster
            expiredSlide.points.stream().filter((o) -> (o.cluster != null)).forEachOrdered((o) -> {
                o.cluster.members.remove(o);
            });
            expiredSlide.points.stream().filter((o) -> (o.cluster == null)).forEachOrdered((o) -> {
                PDList.remove(o);
            });
            //removed expired slide
            window.slides.remove(expiredSlide);
            triggerList.remove(expiredSlide.id);

        }

//        System.out.println("Removing expired slide time = " + (Utils.getCPUTime() - start));
        start = Utils.getCPUTime();

        if (cur_slide_id == Constants.W / Constants.slide) {
            findNeighborFirstWindow();
//                findNeighborNaive();
        } else if (cur_slide_id > Constants.W / Constants.slide) {
            //update neighbor info / cluster for new data points
            Slide newestSlide = window.getNewestSlide();

            use_ub = 0;
            added_to_cluster = 0;

            double foundByUB = 0;
            double foundByLBExact = 0;

            for (ULObject o : newestSlide.points) {
//            newestSlide.points.forEach((ULObject o) -> {
                checked_by_clusters.clear();
                //find cluster using upper bound
//                for (MicroCluster c : microClusters) {
//                    if (UB_distance(o, c.center) <= Constants.R / 2) {
//                        //add o to c
//                        c.members.add(o);
//                        o.cluster = c;
//                        foundByUB +=1;
//                        break;
//                    }
//                    else if(LB_distance(o, c.center) > Constants.R*3.0/2){
//                        checked_by_clusters.addAll(c.members);
//                    }
//                    
//                }

                //find cluster using lb + exact
                if (o.cluster == null) {
                    for (MicroCluster c : microClusters) {
                        if (LB_distance(o, c.center) <= Constants.R / 2) {
//                            System.out.println("R/2 = "+ Constants.R/2);
//                            System.out.println("LB distance = "+LB_distance(o, c.center));
//                            System.out.println("Exact distance = "+ EU_ditance(o, c.center));
                            double exact_distance = EU_ditance(o, c.center);
                            if (exact_distance <= Constants.R / 2) {
                                //add o to c
                                c.members.add(o);
                                o.cluster = c;
                                foundByLBExact += 1;
                                break;

                            } else if (exact_distance > Constants.R * 3 / 2) {
                                checked_by_clusters.addAll(c.members);
                            }
                        }
                    }
                }

                if (o.cluster == null) {
//                    System.out.println("Checked by clusters = "+ checked_by_clusters.size());
                    findNeighbor(o);
                    count += 1;
                } else {
                    added_to_cluster += 1;
                }

            }

            System.out.println("Found by UB = " + foundByUB);
            System.out.println("Found by LB+Exact =" + foundByLBExact);
//            });
        }

        System.out.println("Call find neighbors = " + count * 1.0 / data.size());

//        System.out.println("Added to cluster = " + added_to_cluster);
//        System.out.println("Used upper bound = " + use_ub);
//
//        System.out.println("Finding  neighbor time = " + (Utils.getCPUTime() - start));
//        start = Utils.getCPUTime();
        //update neighbor info of data points in trigger list
        if (cur_slide_id > Constants.W / Constants.slide) {
            if (expiredSlide != null && triggerList.containsKey(expiredSlide.id)) {
                ArrayList<ULObject> triggeredData = triggerList.get(expiredSlide.id);
                triggeredData.stream().filter((o) -> (!o.isSafe)).filter((o) -> (o.neighborCount.containsKey(expiredSlide.id))).map((o) -> {
                    o.numPreEvidence -= o.neighborCount.get(expiredSlide.id);
                    return o;
                }).map((o) -> {
                    o.neighborCount.remove(expiredSlide.id);
                    return o;
                }).filter((o) -> (o.numPreEvidence + o.numSucEvidence < Constants.k)).forEachOrdered((o) -> {
                    reprob(o);
                }); //update neighbor count
            }
        }
//        System.out.println("Re Probinng  neighbor time = " + (Utils.getCPUTime() - start));
        //process dispersed clusters

//        start = Utils.getCPUTime();
        if (cur_slide_id > Constants.W / Constants.slide) {
            for (int i = microClusters.size() - 1; i >= 0; i--) {
                MicroCluster c = microClusters.get(i);
                if (c.members.size() < Constants.k + 1) {
                    //disperse cluster
                    c.members.stream().map((o) -> {
                        o.cluster = null;
                        return o;
                    }).forEachOrdered((o) -> {
                        if (!o.isSafe && o.numPreEvidence + o.numPreEvidence < Constants.k) {
                            reprob(o);
                        }
                    });
                }
                microClusters.remove(i);
            }
        }

//        System.out.println("Dispersed cluster processing time = " + (Utils.getCPUTime() - start));
        //run outlier detection
        PDList.stream().filter((o) -> (o.numPreEvidence + o.numSucEvidence < Constants.k)).forEachOrdered((o) -> {
            outlierList.add(o);
        });
        return outlierList;

    }

}
