package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.PopulationFunction;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author remco@cs.waikato.ac.nz
 * @author Alexei Drummond
 */
@Description("Simulate coalescent trees given a population function.")
public class CoalescentSimulatorBasic extends beast.core.Runnable {
    
    public Input<PopulationFunction> populationFunctionInput
            = new Input<PopulationFunction>(
                    "populationFunction",
                    "the population function",
                    Input.Validate.REQUIRED);
    
    public Input<Integer> replicatesInput = new Input<Integer>(
            "replicates",
            "the number of replicates", Input.Validate.REQUIRED);
    
    public Input<String> outputFileNameInput = new Input<String>(
            "outputFileName",
            "If provided, simulated trees are written to this file rather "
                    + "than to standard out.");

    public Input<TraitSet> timeTraitInput = new Input<TraitSet>(
            "trait",
            "trait information for initializing traits (like node dates) in the tree",
            Input.Validate.REQUIRED);
    
    public Input<Double> maxHeightInput = new Input<Double>(
            "maxHeight",
            "Used to optionally constrain the height of the generated tree "
                    + "to be below the start of the epidemic.  Trees will "
                    + "be generated until the height is found to be smaller "
                    + "than this value.",
            Double.POSITIVE_INFINITY);

    PopulationFunction populationFunction;
    PopulationFunction.Abstract fullPopFunc;

    TraitSet timeTraitSet;
    int replicates;
    String outputFileName;


    @Override
    public void initAndValidate() throws Exception {

        replicates = replicatesInput.get();
        populationFunction = populationFunctionInput.get();
        outputFileName = outputFileNameInput.get();
        timeTraitSet = timeTraitInput.get();

        if (!timeTraitSet.isDateTrait())
            throw new IllegalArgumentException("Trait set must be a date "
                    + "(forward or backward) trait set.");
        
        for (int i=0; i<timeTraitSet.taxaInput.get().asStringList().size(); i++) {
            if (Double.isInfinite(populationFunction.getPopSize(timeTraitSet.getValue(i))))
                throw new IllegalStateException("Intensity is non-finite at "
                        + "one or more sample times.");
        }
        
        
        if (populationFunction instanceof PopulationFunction.Abstract)
            fullPopFunc = (PopulationFunction.Abstract)populationFunction;
    }

    @Override
    public void run() throws Exception {

        PrintStream pstream = null;

        RandomTree tree = new RandomTree();
        for (int i = 0; i < replicates; i++) {
            do {
                if (fullPopFunc != null)
                    fullPopFunc.prepare();
                tree.initByName(
                        "taxonset", timeTraitSet.taxaInput.get(),
                        "trait", timeTraitSet,
                        "populationModel", populationFunction);
            } while (!(tree.getRoot().getHeight()<maxHeightInput.get()));
                    
            // Write output to stdout or file
            if (pstream == null) {
                if (outputFileName == null)
                    pstream = System.out;
                else
                    pstream = new PrintStream(outputFileName);
            }
            
            pstream.print(tree.toString() + ";\n");
        }
        
        if (pstream != null) {
            pstream.flush();
            pstream.close();
        }
    }
    
    
    /**
     * Main method for debugging.
     * 
     * @param args 
     * @throws java.lang.Exception 
     */
    public static void main(String[] args) throws Exception {
        
        int nLeaves = 2;
        int nTrees = 100000;
        
        ConstantPopulation popFun = new ConstantPopulation();
        popFun.initByName("popSize", new RealParameter("1.0"));
        
        List<Taxon> taxonList = new ArrayList<Taxon>();
        for (int i=0; i<nLeaves; i++) {
            String taxonName = "t" + String.valueOf(i+1);
            taxonList.add(new Taxon(taxonName));
        }
        TaxonSet taxonSet = new TaxonSet(taxonList);
        
        StringBuilder traitListBuilder = new StringBuilder();
        for (int i=0; i<taxonList.size(); i++) {
            if (i>0)
                traitListBuilder.append(",");
            
            traitListBuilder.append(taxonList.get(i).getID()).append("=0.0");
        }
        
        TraitSet traitSet = new TraitSet();
        traitSet.initByName(
                "traitname", "date-backward",
                "value", traitListBuilder.toString(),
                "taxa", taxonSet);
        
        CoalescentSimulatorBasic coalSim = new CoalescentSimulatorBasic();
        coalSim.initByName(
                "populationFunction", popFun,
                "replicates", nTrees,
                "outputFileName", "trees.txt",
                "trait", traitSet);
        
        coalSim.run();
    }
}

