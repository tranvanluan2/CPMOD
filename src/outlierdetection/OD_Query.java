/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

/**
 *
 * @author luan
 */
public class OD_Query {

    public double R;
    public int k;
    
    public int W;
    public int S;
//        public int id;

    public OD_Query(double R_, int k_, int W_, int S_) {
        this.R = R_;
        this.k = k_;
        this.W = W_;
        this.S = S_;
    }
}
