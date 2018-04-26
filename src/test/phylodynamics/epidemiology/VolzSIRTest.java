package test.phylodynamics.epidemiology;

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.Coalescent;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.phylodynamics.epidemiology.DeterministicCoalescentSIR;
import beast.phylodynamics.epidemiology.Volz2009TreeDistribution;
import test.beast.BEASTTestCase;

/**
 * @author Alexei Drummond
 */
public class VolzSIRTest extends BEASTTestCase {
    String[] trees = new String[]{"((A:1.0,B:1.0):1.0,C:2.0);"}; //, "" more trees ?
    Alignment data;
    final double n_S0 = 100;
    final double beta = 0.02;
    final double gamma = 0.5;
    final double origin = 5;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        data = getFourTaxaNoData();
    }

    public void testVolzSIR() throws Exception {
        // *********** 3 taxon **********
        Tree tree = getTree(data, trees[0]);
        TreeIntervals treeIntervals = new TreeIntervals();
        treeIntervals.initByName("tree", tree);

        DeterministicCoalescentSIR detCoalSIR = new DeterministicCoalescentSIR();
        detCoalSIR.initByName(
                "n_S0", Double.toString(n_S0),
                "beta", Double.toString(beta),
                "gamma", Double.toString(gamma),
                "origin", Double.toString(origin),
                "integrationStepCount", 1000);

        // TODO: failed because detCoalSIR is not inherited from PopulationFunction,
        Coalescent coal = new Coalescent();
        coal.initByName("treeIntervals", treeIntervals, "populationModel", detCoalSIR);
//        SIRPopulationFunction popFunction = new SIRPopulationFunction();
//        popFunction.initByName("volzSIR", detCoalSIR);

        double logL = coal.calculateLogP();
        // TODO: calculateLogLikelihood not implemented
//        double logL = detCoalSIR.calculateLogLikelihood(treeIntervals, popFunction);
        System.out.println("Coalescent(VolzSIR).logL=" + logL);

        Volz2009TreeDistribution volzCoalescent = new Volz2009TreeDistribution();
        volzCoalescent.initByName(
                "S0", Double.toString(n_S0),
                "beta", Double.toString(beta),
                "delta", Double.toString(gamma),
                "T", Double.toString(origin),
                "treeIntervals", treeIntervals);

        double logL2 = volzCoalescent.calculateLogP();

        assertEquals(logL, logL2, PRECISION);
    }

    public void testExponentialGrowth() throws Exception {

    }

}
