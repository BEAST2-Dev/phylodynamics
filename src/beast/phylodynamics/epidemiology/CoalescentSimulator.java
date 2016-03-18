package beast.phylodynamics.epidemiology;

/**
 * @author Alexei Drummond
 */

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
//import beast.evolution.tree.TreeTraceAnalysis;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author remco@cs.waikato.ac.nz
 */
@Description("Performs random tree generation for a given site model. ")
public class CoalescentSimulator extends beast.core.Runnable {
    public Input<PopulationFunction> populationFunctionInput
            = new Input<PopulationFunction>("populationFunction", "the population function", Input.Validate.REQUIRED);
    public Input<Integer> ntaxaInput = new Input<Integer>("ntaxa", "the number of taxa", Input.Validate.REQUIRED);
    public Input<Integer> replicatesInput = new Input<Integer>("replicates", "the number of replicates", Input.Validate.REQUIRED);
    public Input<String> outputFileNameInput = new Input<String>(
            "outputFileName",
            "If provided, simulated trees are written to this file rather "
                    + "than to standard out.");

    public Input<TaxonSet> m_taxonset = new Input<TaxonSet>("taxonset","set of taxa to initialise tree with specified by a taxonset", Input.Validate.REQUIRED);

    public Input<List<TraitSet>> m_traitList = new Input<List<TraitSet>>("trait",
            "trait information for initializing traits (like node dates) in the tree",
            new ArrayList<TraitSet>(), Input.Validate.REQUIRED);

    /**
     * nr of categories in site model *
     */
    int ntaxa;
    /**
     * nr of states in site model *
     */
    int replicates;

    PopulationFunction populationFunction;

    TaxonSet taxa = new TaxonSet();

    List<TraitSet> timeTraitSet;

    /**
     * name of output file *
     */
    String outputFileName;


    @Override
    public void initAndValidate() {

        replicates = replicatesInput.get();
        ntaxa = ntaxaInput.get();
        populationFunction = populationFunctionInput.get();
        outputFileName = outputFileNameInput.get();
        timeTraitSet = m_traitList.get();

        List<Taxon> taxaList = new ArrayList<Taxon>();
        for (int i = 1; i <= ntaxa; i++) {
            taxaList.add(new Taxon(i + ""));
        }
        taxa.initByName("taxon", taxaList);
    }

    @Override
    public void run() throws Exception {

        // Write output to stdout or file
        PrintStream pstream;
        if (outputFileName == null)
            pstream = System.out;
        else
            pstream = new PrintStream(outputFileName);

        List<Tree> treeAnalysis = new ArrayList<Tree>();

        for (int i = 0; i < replicates; i++) {
            RandomTree tree = new RandomTree();
            tree.initByName("taxonset", taxa, "trait", timeTraitSet, "populationModel", populationFunction);

            treeAnalysis.add(tree);
            pstream.print(tree.toString() + "\n");
        }

        // Use TreeTraceAnalysis to report topologies
        //TreeTraceAnalysis tta = new TreeTraceAnalysis(treeAnalysis, 0.0, 95.0);
        //tta.report(pstream);

        if (outputFileName != null) {
            pstream.flush();
            pstream.close();
        }
    }
} // class CoalescentSimulator

