/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package test.phylodynamics.epidemiology;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.phylodynamics.epidemiology.StochasticCoalescent;
import beast.phylodynamics.epidemiology.StochasticSIR;
import beast.phylodynamics.epidemiology.StochasticSIRPopulationFunction;
import beast.util.TreeParser;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class StochasticCoalescentTest {
    
    public StochasticCoalescentTest() { }

    /**
     * Test that the tree density conditioned on a single trajectory is
     * calculated correctly.
     * 
     * @throws java.lang.Exception because BEASTObject.initByName throws
     * this.
     */
    @Test
    public void trajectoryConditionedLikelihood() throws Exception {

        double origin = 12.7808530307;
        double S0 = 999;
        double beta = 0.00075;
        double gamma = 0.4;
        int latticeSize = 1000;
        
        // Set up trajectory object
        StochasticSIR stochasticSIR = new StochasticSIR();
        
        stochasticSIR.initByName(
                "n_S0", new RealParameter(String.valueOf(S0)),
                "beta", new RealParameter(String.valueOf(beta)),
                "gamma", new RealParameter(String.valueOf(gamma)),
                "origin", new RealParameter(String.valueOf(origin)),
                "numSamplesFromTrajectory", latticeSize,
                "integrationStepCount", latticeSize);

        /*
        Generate trajectories until stochasticSIR is no longer "dirty". 
        This is a hack, becaues we don't actually care about this simulated
        trajectory - we have our own.  However, we require the dirty flag
        to be set to false, which doesn't happen unless a "good" trajectory
        is generated.
        */
        while(!stochasticSIR.simulateTrajectory());
        
        // Load in pre-generated trajectory from file
        List<Double> trajTime = new ArrayList<Double>();
        List<Double> trajS = new ArrayList<Double>();
        List<Double> trajI = new ArrayList<Double>();
        
        InputStream in = this.getClass().getResourceAsStream("TestTrajectory.txt");
        BufferedReader breader = new BufferedReader(new InputStreamReader(in));
        String line;
        boolean isFirst = true;
        while ((line = breader.readLine()) != null) {
            if (isFirst) {
                // Skip header
                isFirst = false;
                continue;
            }
            
            String [] strVals = line.trim().split(" ");
            if (strVals.length==3) {
                trajTime.add(Double.valueOf(strVals[0]));
                trajS.add(Double.valueOf(strVals[1]));
                trajI.add(Double.valueOf(strVals[2]));
            }
        }
        breader.close();
        
        // Lay trajectory values onto regularly spaced lattice
        int latticeTime = 0;
        double latticeDt = origin/(latticeSize-1);
        int latticeIdx = 0;
        
        stochasticSIR.NStraj.clear();
        stochasticSIR.NItraj.clear();
        
        for (int eventIdx=0; eventIdx<trajTime.size(); eventIdx++) {
            
            double tNext;
            if (eventIdx<trajTime.size()-1)
                tNext = trajTime.get(eventIdx+1);
            else
                tNext = origin;
            
            while (latticeIdx*latticeDt<tNext) {
                stochasticSIR.NStraj.add(trajS.get(eventIdx));
                stochasticSIR.NItraj.add(trajI.get(eventIdx));
                latticeIdx += 1;
            }
            
        }
        
        // Make sure trajectory at last lattice point is recorded
        if (stochasticSIR.NStraj.size()<latticeSize) {
            stochasticSIR.NStraj.add(trajS.get(trajS.size()-1));
            stochasticSIR.NItraj.add(trajI.get(trajI.size()-1));
        }
        
        // DEBUG: Write trajectory lattice to disk
//        PrintStream outf = new PrintStream("traj.txt");
//        outf.println("t S I");
//        for (int i=0; i<latticeSize; i++) {
//            outf.println(i*latticeDt + " " + stochasticSIR.NStraj.get(i) + " " + stochasticSIR.NItraj.get(i));
//        }
//        outf.close();
        
        // Calculate "effective" population sizes:
        stochasticSIR.effectivePopSizeTraj.clear();
        
        for (int i=0; i<latticeSize; i++) {
            double S = stochasticSIR.NStraj.get(i);
            double I = stochasticSIR.NItraj.get(i);
            
            stochasticSIR.effectivePopSizeTraj.add((I-1)/(2.0*beta*S));
        }
        
        // Switch effective pop size to reverse time
        Collections.reverse(stochasticSIR.effectivePopSizeTraj);

        // Estimate intensity on integration lattice:
        stochasticSIR.intensityTraj.clear();

        double intensity = 0.0;
        stochasticSIR.intensityTraj.add(intensity);

        for (Double popSize : stochasticSIR.effectivePopSizeTraj) {
            intensity += latticeDt / popSize;
            stochasticSIR.intensityTraj.add(intensity);
        }

        // Start of integral is 0.5*dt from end of forward-time integration.
        stochasticSIR.tIntensityTrajStart = -0.5 * latticeDt;
        
        // Set up tree
        String newickStr = "(((79:0.642827446217,87:1.56573879247):0.174586994043,(22:2.2990314372,(56:2.95219148906,((85:1.86532867775,(7:0.463894038781,(47:2.31882898943,(43:0.863163753429,(41:0.125110395388,36:0.642400616991):0.0449054318605):0.7072548983):1.41719614874):0.647968801804):0.74692121537,(21:0.329699281024,(((95:0.209230044325,65:1.02709362659):1.73124264848,((27:0.551768147649,5:0.157775910655):2.44821236615,((82:1.20237560354,48:3.03421141731):1.1621193946,(51:0.107013874421,9:1.36130638469):6.66312891406):0.113967142007):1.49812797593):0.0963213204681,(8:0.456361390405,(80:0.341589758911,(88:0.409882092896,((92:0.509324611178,(58:1.25046032829,83:0.221875292397):0.12040350037):4.16976137118,(63:0.705682501185,(29:1.23886194456,(90:1.23063843486,42:2.02595751112):0.510722272655):0.400296332537):2.30902513584):0.203813512637):0.358155731828):0.419950767736):3.23909167044):0.274519724595):0.283461892955):0.289213482822):0.364262389046):0.501668803208):1.00139330805,(50:2.64586140577,(6:2.1245900793,(((96:0.389740462205,10:1.24965068456):0.752917336836,((55:4.78101153052,17:5.85888008943):3.57188264344,(((77:3.1958770881,91:0.321318157708):1.17360302095,(52:1.61631375224,((23:1.02441580121,35:0.557274287256):0.346294678688,(30:1.57804157244,(11:0.116605937481,98:0.288425921491):0.546000064604):1.4142945411):1.71086142359):0.800268333908):2.42619450434,((100:0.363618051948,(18:4.26784371269,(32:0.135617350568,(94:1.18531998245,76:1.66330444203):1.50372564217):0.681000607715):1.33956017126):2.48809274104,(37:0.275643035649,((70:1.29592636384,(97:1.40208031833,(25:1.59478677769,1:5.94136906579):0.0903138317307):0.104931213787):0.116015269773,(4:3.7278908785,((33:2.35689629239,60:0.0185789162998):1.67579414906,(14:3.10968745891,(24:1.33934014401,39:0.366322121588):1.01917004066):0.131875272431):0.000245438325695):0.0356492297219):1.75607035758):0.184952384951):0.801189642948):0.067110682324):0.467621494323):0.89621212524,(((93:2.41623901896,(86:0.657870063822,75:0.74103527635):5.00796218632):0.00619596972062,(45:2.6021275173,(3:0.0831178528513,20:2.55657590577):0.527074975774):1.94612822902):2.40980124466,((74:1.88373467917,((31:0.820128559822,(49:0.773253848525,99:0.153300427073):0.143847475813):0.66799866422,((19:0.257152601127,(64:2.78877239902,59:0.118112118434):0.242501954033):0.150613951396,((62:0.172136897933,(73:0.785157966353,44:1.56998679617):2.77376987512):2.02263798582,((78:1.28104215284,(66:0.480074650411,61:1.004735098):0.897849464716):1.77914662441,(72:0.650335241107,(28:0.644896407049,46:1.97504344831):2.2523487189):0.487313396696):0.313413755117):0.348426554545):0.0814293534687):0.496395953197):0.747329556311,((67:0.585945255164,40:0.209636710833):0.289687605067,(69:0.526019884768,((53:0.548354635462,(57:1.53527087874,(34:0.801752408612,(12:1.2182266603,(81:1.10143886873,26:0.73649039298):0.843544121588):1.65040665807):0.0185782900061):1.28075539663):1.06907735375,((68:0.749231041113,(16:0.552374881678,(38:0.536360899103,(84:2.13535735464,71:2.87728905599):0.203647605566):0.0942290019581):0.577552249027):0.907395294241,(2:1.79371173528,(89:0.394869633041,(15:3.17634334689,(13:2.56757595841,54:2.81147606244):0.062178166339):0.198354387682):0.0201536146416):1.25099023532):0.908408820621):0.173944559366):2.41148973972):0.41096364666):0.709591312145):1.07270721056):0.872885640466):0.420419827939):0.397480218169):0.295470991713;";
        TreeParser tree = new TreeParser();
        tree.initByName(
                "adjustTipHeights", false,
                "newick", newickStr);
        TreeIntervals treeIntervals = new TreeIntervals(tree);
        
        // Set up 
        StochasticCoalescent stochasticCoalescent = new StochasticCoalescent();
        StochasticSIRPopulationFunction popFunction =
                new StochasticSIRPopulationFunction();
        popFunction.initByName("stochasticSIR", stochasticSIR);

        // Calculate probability denstiy of tree given trajectory
        double logP = stochasticCoalescent.calculateLogLikelihood(
                treeIntervals, popFunction);
        
        double logPtruth = -1763.206425;  // from likelihoodEstimator.R
        
        System.out.println("logP: " + logP);
        System.out.println("True logP:" + logPtruth);
        
        // Success if <1% error
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-2);
    }
}
