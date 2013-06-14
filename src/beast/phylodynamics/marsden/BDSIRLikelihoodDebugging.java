/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package beast.phylodynamics.marsden;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.phylodynamics.epidemiology.HybridSEIREpidemic;
import beast.util.TreeParser;
import java.io.PrintStream;

/**
 * Code I'm using to try to understand what BDSIR is doing.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BDSIRLikelihoodDebugging {
    
    public static void main (String [] args) throws Exception {
        
        // Initialise tree
        String newick = "((20:1.7327013309920987,(2:1.5699517642369987,9:1.5699517642369987):0.1627495667551):0.3906420437073396,(((((12:1.4606660366895723,((22:0.7364772016601255,21:0.7364772016601255):0.6701138844558105,(15:1.3665082235172927,((3:1.020945860930602,(8:0.8653300121493379,10:0.8653300121493379):0.155615848781264):0.22753196406186138,(7:1.1290942502849537,19:1.1290942502849537):0.1193835747075096):0.11803039852482944):0.040082862598643354):0.054074950573636205):0.0074844678674776954,((17:0.012266430667138728,23:0.012266430667138728):1.3007893648454143,(24:0.9390090648133209,5:0.9390090648133209):0.3740467306992321):0.15509470904449696):0.038425441676050665,((18:0.5287382111016794,16:0.5287382111016794):0.6219873839148018,14:1.1507255950164812):0.35585035121661945):0.17026942978592519,11:1.6768453760190258):0.3913503628622723,((13:0.7724098713155667,(25:0.07463738743032255,4:0.07463738743032255):0.6977724838852442):0.9249671690996732,(6:1.2765864559055629,1:1.2765864559055629):0.4207905845096771):0.3708186984660582):0.05514763581814014);";
        Tree tree = new TreeParser(newick, false);
        TreeIntervals treeIntervals = new TreeIntervals();
        treeIntervals.init(tree);
        
        // Generate ensemble of SIR trajectories
        HybridSEIREpidemic sir = new HybridSEIREpidemic();
        sir.initByName(
                "tree", tree,
                "birth", new RealParameter("0.005"),
                "death", new RealParameter("2.0"),
                "sampling", "0",
                "origin", new RealParameter("2.5"),
                "S0", new RealParameter("1000.0"),
                "Nt", "10001",
                "Nsamples", "101",
                "simulationType", "SAL", "isDeterministic", false);
        
        int nTraj = 1000;
        double [] meanS = new double[10001];
        double [] varS = new double[10001];
        
        // Simulate trajectories:
        for (int traj=0; traj<nTraj; traj++) {
            sir.initTraj(0.005, 0., new Double[]{2.0}, new Double[]{0.}, 2.5, tree.getLeafNodeCount(), 100, 10001, sir.times);

            double S = 1000.0;
            meanS[0] += S;
            varS[0] += S*S;
            for (int i=0; i<sir.dS.length; i++) {
                S -= sir.dS[i];
                meanS[i+1] += S;
                varS[i+1] += S*S;
            }
        }
        
        // Normalize results:
        for (int i=0; i<sir.dS.length+1; i++) {
            meanS[i] /= nTraj;
            varS[i] = varS[i]/nTraj - meanS[i]*meanS[i];
        }
        
        // Write to disk:
        double dt = 2.5/(sir.dS.length);
        PrintStream outFile = new PrintStream("BDSIRLikelihoodDebugging.txt");
        outFile.println("t meanS varS");
        for (int i=0; i<sir.dS.length+1; i++) {
            outFile.format("%g %g %g\n", dt*i, meanS[i], varS[i]);
        }
        outFile.close();
    }
}
