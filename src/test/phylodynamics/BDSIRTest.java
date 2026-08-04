package test.phylodynamics;

import beast.base.core.Description;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;
import phylodynamics.BDSIR;
import beast.base.util.Randomizer;

/**
 * @author Shubhangi
 */

@Description("Test for comparison during BEAST3 migration. logPs computed from BDSIR v1.5 based on bdsky 1.5 ")
public class BDSIRTest {

    public BDSIRTest() { }
	/**
	 * Basic test for BDSIR
	 * Reference from phylodynamics version 1.5.0
	 */
    @Test
    public void testBDSIRLikelihood() throws Exception {

        Randomizer.setSeed(42);
        String newickStr = "(((22:2.2990314372014935,(56:2.9521914890587606,((21:0.32969928102424584,(((95:0.20923004432474013,65:1.0270936265943202):1.731242648480694,((27:0.5517681476494722,5:0.15777591065467167):2.448212366152161,((82:1.2023756035356818,48:3.034211417314549):1.1621193946048294,(51:0.10701387442070143,9:1.3613063846928863):6.663128914063232):0.11396714200741442):1.4981279759258852):0.09632132046807573,(8:0.4563613904047621,((((63:0.7056825011848105,((90:1.2306384348590065,42:2.025957511118051):0.510722272654716,29:1.2388619445646505):0.4002963325371933):2.309025135836837,((58:1.2504603282882396,83:0.22187529239694292):0.12040350036963865,92:0.5093246111777496):4.169761371182936):0.20381351263653613,88:0.40988209289627164):0.35815573182832416,80:0.34158975891102905):0.41995076773605167):3.239091670436549):0.2745197245947222):0.2834618929550565,((7:0.4638940387809063,(47:2.318828989431455,(43:0.8631637534291592,(41:0.12511039538846092,36:0.642400616991365):0.04490543186047358):0.7072548982999889):1.417196148739511):0.6479688018042475,85:1.8653286777482379):0.7469212153703846):0.28921348282221304):0.3642623890464143):0.5016688032075132,(79:0.6428274462169905,87:1.5657387924674513):0.1745869940428697):1.0013933080542134,((((((((100:0.36361805194798436,((32:0.13561735056803315,(94:1.1853199824473748,76:1.6633044420319116):1.5037256421720855):0.6810006077145765,18:4.267843712685522):1.3395601712634475):2.488092741044201,(((((25:1.594786777694054,1:5.941369065789077):0.09031383173071283,97:1.4020803183322919):0.10493121378674086,70:1.2959263638436456):0.11601526977289112,(((33:2.3568962923913546,60:0.018578916299765602):1.6757941490565083,(14:3.1096874589082715,(24:1.3393401440122652,39:0.36632212158815936):1.0191700406631137):0.13187527243101016):0.0002454383256953463,4:3.72789087849912):0.03564922972192086):1.7560703575750258,37:0.2756430356489643):0.18495238495119537):0.8011896429483154,((77:3.1958770881049485,91:0.32131815770842564):1.1736030209546113,(((23:1.0244158012059774,35:0.5572742872564813):0.34629467868753494,((11:0.11660593748062809,98:0.2884259214909246):0.5460000646040939,30:1.5780415724425616):1.4142945411023877):1.7108614235922541,52:1.6163137522421795):0.8002683339078454):2.4261945043445476):0.06711068232398532,(55:4.781011530516491,17:5.858880089428844):3.5718826434382662):0.46762149432250366,(96:0.38974046220477154,10:1.2496506845631044):0.7529173368363931):0.8962121252404045,(((74:1.8837346791718197,(((19:0.2571526011271228,(64:2.7887723990200106,59:0.11811211843360692):0.2425019540326998):0.15061395139649658,((62:0.1721368979327691,(73:0.7851579663529034,44:1.5699867961689673):2.7737698751218858):2.0226379858191086,((78:1.2810421528437788,(66:0.4800746504106126,61:1.0047350980045788):0.8978494647157946):1.779146624408618,((28:0.644896407049,46:1.975043448314823):2.25234871889561,72:0.6503352411065988):0.4873133966963845):0.3134137551166951):0.34842655454507465):0.08142935346869962,(31:0.8201285598215611,(49:0.7732538485247655,99:0.15330042707285685):0.14384747581345358):0.6679986642198932):0.4963959531967346):0.7473295563109597,((((53:0.5483546354617461,(57:1.5352708787432974,(34:0.8017524086121206,(12:1.2182266602962883,(81:1.101438868727076,26:0.7364903929803983):0.8435441215878487):1.6504066580694872):0.018578290006075804):1.2807553966307017):1.0690773537474865,(((16:0.5523748816782845,(38:0.5363608991026343,(84:2.1353573546402913,71:2.8772890559882):0.20364760556639894):0.09422900195806072):0.5775522490272174,68:0.7492310411126208):0.9073952942408914,((89:0.3948696330409902,((13:2.5675759584117337,54:2.8114760624370696):0.062178166339046825,15:3.1763433468949263):0.19835438768216918):0.020153614641559514,2:1.7937117352805494):1.2509902353242888):0.9084088206213767):0.1739445593663662,69:0.5260198847676936):2.4114897397231223,(67:0.5859452551639297,40:0.20963671083264845):0.28968760506734803):0.41096364665959584):0.7095913121446329,((93:2.416239018964405,(86:0.6578700638218109,75:0.7410352763503614):5.007962186321123):0.006195969720623751,((3:0.08311785285125595,20:2.556575905767927):0.5270749757736688,45:2.602127517303609):1.946128229019207):2.4098012446604966):1.0727072105608606):0.8728856404660741,6:2.124590079296869):0.4204198279390081,50:2.6458614057743235):0.39748021816855017):0.29547099171304114;";
        TreeParser tree = new TreeParser();
        tree.initByName("adjustTipHeights", false, "newick", newickStr);

        // Set up trajectory object
        BDSIR bdsir = new BDSIR();

        bdsir.initByName(
                "S0", new RealParameter(String.valueOf(1000.)),
                "tree", tree,
                "reproductiveNumber", new RealParameter(String.valueOf(2.5)),
                "becomeUninfectiousRate", new RealParameter(String.valueOf(0.2)),
                "samplingProportion", new RealParameter(String.valueOf(0.10)),
                "dS", new RealParameter("100. 200."),
                "dR", new RealParameter("10. 20."),
                "origin", new RealParameter(String.valueOf(15.)));

        double logP = bdsir.calculateTreeLogLikelihood(tree);
        double logPtruth = -451.3896063922145;
        System.out.println("Scenario 1: logP=" + logP + " truth=" + logPtruth);
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-4);


        bdsir.initByName(
                "S0", new RealParameter(String.valueOf(10000.)),
                "tree", tree,
                "reproductiveNumber", new RealParameter(String.valueOf(2.5)),
                "becomeUninfectiousRate", new RealParameter(String.valueOf(0.2)),
                "samplingProportion", new RealParameter(String.valueOf(0.10)),
                "dS", new RealParameter("1000. 1000."),
                "dR", new RealParameter("100. 200."),
                "origin", new RealParameter(String.valueOf(15.)));

        logP = bdsir.calculateTreeLogLikelihood(tree);
		// logP computed from BDSIR from phylodynamics version 1.5.0 based on bdsky version 1.5.1 for comparison during BEAST3 migration
        logPtruth = -451.2857974679011;
        System.out.println("Scenario 2: logP=" + logP + " truth=" + logPtruth);
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-4);


        newickStr = "(((((1:23.201674279943926,3:21.201674279943926):2.7431760554195783,14:12.944850335363505):1.2575389333037634,((4:3.3051116149767834,7:0.3051116149767834):18.02597361151068,20:5.331085226487463):2.871304042179805):0.6103254377900278,((((5:20.875764555662627,12:13.875764555662627):0.4167891392779133,17:9.29255369494054):1.63903996381406,(10:9.766039351017923,16:3.7660393510179233):8.165554307736677):0.7197286663179412,(6:15.933403099522199,15:6.933403099522199):6.717919225550343):0.1613923813847542):0.3310475410146694,((((2:20.97070734376335,19:3.970707343763351):1.3074622585163276,(8:9.41339312649908,13:4.413393126499081):6.864776475780598):0.8823641166101801,18:7.160533718889859):1.0497844441524649,(9:4.853545654572189,11:2.853545654572189):12.356772508470135):2.9334440844296417):0.0;";
        tree = new TreeParser();
        tree.initByName("adjustTipHeights", false, "newick", newickStr);
        bdsir.initByName(
                "S0", new RealParameter(String.valueOf(1000.)),
                "tree", tree,
                "reproductiveNumber", new RealParameter(String.valueOf(2.5)),
                "becomeUninfectiousRate", new RealParameter(String.valueOf(0.2)),
                "samplingProportion", new RealParameter(String.valueOf(0.10)),
                "dS", new RealParameter("100. 200."),
                "dR", new RealParameter("10. 20."),
                "origin", new RealParameter(String.valueOf(30.)));

        logP = bdsir.calculateTreeLogLikelihood(tree);
		// logP computed from BDSIR from phylodynamics version 1.5.0 based on bdsky version 1.5.1 for comparison during BEAST3 migration
        logPtruth = -157.40908231942964;
        System.out.println("Scenario 3: logP=" + logP + " truth=" + logPtruth);
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-4);


        bdsir.initByName(
                "S0", new RealParameter(String.valueOf(1000.)),
                "tree", tree,
                "reproductiveNumber", new RealParameter(String.valueOf(5)),
                "becomeUninfectiousRate", new RealParameter(String.valueOf(0.2)),
                "samplingProportion", new RealParameter(String.valueOf(0.01)),
                "dS", new RealParameter("100. 200."),
                "dR", new RealParameter("10. 20."),
                "origin", new RealParameter(String.valueOf(30.)));

        logP = bdsir.calculateTreeLogLikelihood(tree);
		// logP computed from BDSIR from phylodynamics version 1.5.0 based on bdsky version 1.5.1 for comparison during BEAST3 migration
        logPtruth = -278.5405119276242;
        System.out.println("Scenario 4: logP=" + logP + " truth=" + logPtruth);
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-4);


        bdsir.initByName(
                "S0", new RealParameter(String.valueOf(1000.)),
                "tree", tree,
                "reproductiveNumber", new RealParameter(String.valueOf(1.5)),
                "becomeUninfectiousRate", new RealParameter(String.valueOf(0.2)),
                "samplingProportion", new RealParameter(String.valueOf(0.10)),
                "dS", new RealParameter("100. 200."),
                "dR", new RealParameter("10. 20."),
                "origin", new RealParameter(String.valueOf(50.)));

        logP = bdsir.calculateTreeLogLikelihood(tree);
		// logP computed from BDSIR from phylodynamics version 1.5.0 based on bdsky version 1.5.1 for comparison during BEAST3 migration
        logPtruth = -132.03873262406228;
        System.out.println("Scenario 5: logP=" + logP + " truth=" + logPtruth);
        assert(Math.abs(logP-logPtruth)/Math.abs(0.5*(logP+logPtruth))<1e-4);
    }
}