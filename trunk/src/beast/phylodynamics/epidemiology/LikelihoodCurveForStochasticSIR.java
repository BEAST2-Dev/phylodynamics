package beast.phylodynamics.epidemiology;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import beast.evolution.tree.coalescent.*;

import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author Alex Popinga
 * @author Alexei Drummond
 */
public class LikelihoodCurveForStochasticSIR {

    public static final String TRUE_TREE = "(((22[NEWICKtype='I',reaction='Recovery',time=4.097564540176261]:2.2990314372014935,(56[NEWICKtype='I',reaction='Sampling',time=5.114986981079943]:2.9521914890587606,((21[NEWICKtype='I',reaction='Recovery',time=3.0651701488226974]:0.32969928102424584,(((95[NEWICKtype='I',reaction='Recovery',time=5.046784605666684]:0.20923004432474013,65[NEWICKtype='I',reaction='Recovery',time=5.864648187936264]:1.0270936265943202)[NEWICKtype='I',reaction='Infection',time=4.8375545613419435]:1.731242648480694,((27[NEWICKtype='I',reaction='Recovery',time=7.604420402588768]:0.5517681476494722,5[NEWICKtype='I',reaction='Recovery',time=7.210428165593967]:0.15777591065467167)[NEWICKtype='I',reaction='Infection',time=7.052652254939296]:2.448212366152161,((82[NEWICKtype='I',reaction='Recovery',time=7.08290202893506]:1.2023756035356818,48[NEWICKtype='I',reaction='Recovery',time=8.914737842713928]:3.034211417314549)[NEWICKtype='I',reaction='Infection',time=5.8805264253993785]:1.1621193946048294,(51[NEWICKtype='I',reaction='Sampling',time=11.488549819278482]:0.10701387442070143,9[NEWICKtype='I',reaction='Sampling',time=12.742842329550667]:1.3613063846928863)[NEWICKtype='I',reaction='Infection',time=11.38153594485778]:6.663128914063232)[NEWICKtype='I',reaction='Infection',time=4.718407030794549]:0.11396714200741442)[NEWICKtype='I',reaction='Infection',time=4.604439888787135]:1.4981279759258852)[NEWICKtype='I',reaction='Infection',time=3.1063119128612495]:0.09632132046807573,(8[NEWICKtype='I',reaction='Sampling',time=6.705443653234485]:0.4563613904047621,((((63[NEWICKtype='I',reaction='Recovery',time=10.245709912052282]:0.7056825011848105,((90[NEWICKtype='I',reaction='Recovery',time=11.681684450918388]:1.2306384348590065,42[NEWICKtype='I',reaction='Recovery',time=12.477003527177432]:2.025957511118051)[NEWICKtype='I',reaction='Infection',time=10.451046016059381]:0.510722272654716,29[NEWICKtype='I',reaction='Recovery',time=11.179185687969316]:1.2388619445646505)[NEWICKtype='I',reaction='Infection',time=9.940323743404665]:0.4002963325371933)[NEWICKtype='I',reaction='Infection',time=9.540027410867472]:2.309025135836837,((58[NEWICKtype='I',reaction='Recovery',time=12.77162747487145]:1.2504603282882396,83[NEWICKtype='I',reaction='Recovery',time=11.743042438980153]:0.22187529239694292)[NEWICKtype='I',reaction='Infection',time=11.52116714658321]:0.12040350036963865,92[NEWICKtype='I',reaction='Recovery',time=11.91008825739132]:0.5093246111777496)[NEWICKtype='I',reaction='Infection',time=11.400763646213571]:4.169761371182936)[NEWICKtype='I',reaction='Infection',time=7.231002275030635]:0.20381351263653613,88[NEWICKtype='I',reaction='Recovery',time=7.43707085529037]:0.40988209289627164)[NEWICKtype='I',reaction='Infection',time=7.027188762394099]:0.35815573182832416,80[NEWICKtype='I',reaction='Recovery',time=7.010622789476804]:0.34158975891102905)[NEWICKtype='I',reaction='Infection',time=6.669033030565775]:0.41995076773605167)[NEWICKtype='I',reaction='Infection',time=6.249082262829723]:3.239091670436549)[NEWICKtype='I',reaction='Infection',time=3.009990592393174]:0.2745197245947222)[NEWICKtype='I',reaction='Infection',time=2.7354708677984516]:0.2834618929550565,((7[NEWICKtype='I',reaction='Recovery',time=4.3107930307989335]:0.4638940387809063,(47[NEWICKtype='I',reaction='Recovery',time=7.582924130188993]:2.318828989431455,(43[NEWICKtype='I',reaction='Recovery',time=6.834513792486686]:0.8631637534291592,(41[NEWICKtype='I',reaction='Recovery',time=6.141365866306462]:0.12511039538846092,36[NEWICKtype='I',reaction='Recovery',time=6.658656087909366]:0.642400616991365)[NEWICKtype='I',reaction='Infection',time=6.016255470918001]:0.04490543186047358)[NEWICKtype='I',reaction='Infection',time=5.971350039057527]:0.7072548982999889)[NEWICKtype='I',reaction='Infection',time=5.264095140757538]:1.417196148739511)[NEWICKtype='I',reaction='Infection',time=3.846898992018027]:0.6479688018042475,85[NEWICKtype='I',reaction='Recovery',time=5.0642588679620175]:1.8653286777482379)[NEWICKtype='I',reaction='Infection',time=3.1989301902137797]:0.7469212153703846)[NEWICKtype='I',reaction='Infection',time=2.452008974843395]:0.28921348282221304)[NEWICKtype='I',reaction='Infection',time=2.162795492021182]:0.3642623890464143)[NEWICKtype='I',reaction='Infection',time=1.7985331029747678]:0.5016688032075132,(79[NEWICKtype='I',reaction='Recovery',time=2.114278740027115]:0.6428274462169905,87[NEWICKtype='I',reaction='Recovery',time=3.0371900862775756]:1.5657387924674513)[NEWICKtype='I',reaction='Infection',time=1.4714512938101243]:0.1745869940428697)[NEWICKtype='I',reaction='Infection',time=1.2968642997672546]:1.0013933080542134,((((((((100[NEWICKtype='I',reaction='Sampling',time=7.070101416114068]:0.36361805194798436,((32[NEWICKtype='I',reaction='Recovery',time=8.86266149371214]:0.13561735056803315,(94[NEWICKtype='I',reaction='Recovery',time=11.416089767763568]:1.1853199824473748,76[NEWICKtype='I',reaction='Recovery',time=11.894074227348105]:1.6633044420319116)[NEWICKtype='I',reaction='Infection',time=10.230769785316193]:1.5037256421720855)[NEWICKtype='I',reaction='Infection',time=8.727044143144107]:0.6810006077145765,18[NEWICKtype='I',reaction='Recovery',time=12.313887248115053]:4.267843712685522)[NEWICKtype='I',reaction='Infection',time=8.046043535429531]:1.3395601712634475)[NEWICKtype='I',reaction='Infection',time=6.7064833641660835]:2.488092741044201,(((((25[NEWICKtype='I',reaction='Recovery',time=8.065460458632502]:1.594786777694054,1[NEWICKtype='I',reaction='Recovery',time=12.412042746727526]:5.941369065789077)[NEWICKtype='I',reaction='Infection',time=6.470673680938448]:0.09031383173071283,97[NEWICKtype='I',reaction='Sampling',time=7.782440167540027]:1.4020803183322919)[NEWICKtype='I',reaction='Infection',time=6.380359849207736]:0.10493121378674086,70[NEWICKtype='I',reaction='Recovery',time=7.57135499926464]:1.2959263638436456)[NEWICKtype='I',reaction='Infection',time=6.275428635420995]:0.11601526977289112,(((33[NEWICKtype='I',reaction='Recovery',time=10.227998475143583]:2.3568962923913546,60[NEWICKtype='I',reaction='Recovery',time=7.889681099051994]:0.018578916299765602)[NEWICKtype='I',reaction='Infection',time=7.871102182752228]:1.6757941490565083,(14[NEWICKtype='I',reaction='Recovery',time=9.436870765035001]:3.1096874589082715,(24[NEWICKtype='I',reaction='Recovery',time=8.685693490802109]:1.3393401440122652,39[NEWICKtype='I',reaction='Recovery',time=7.712675468378003]:0.36632212158815936)[NEWICKtype='I',reaction='Infection',time=7.346353346789844]:1.0191700406631137)[NEWICKtype='I',reaction='Infection',time=6.32718330612673]:0.13187527243101016)[NEWICKtype='I',reaction='Infection',time=6.19530803369572]:2.454383256953463E-4,4[NEWICKtype='I',reaction='Sampling',time=9.922953473869145]:3.72789087849912)[NEWICKtype='I',reaction='Infection',time=6.1950625953700245]:0.03564922972192086)[NEWICKtype='I',reaction='Infection',time=6.159413365648104]:1.7560703575750258,37[NEWICKtype='I',reaction='Recovery',time=4.678986043722042]:0.2756430356489643)[NEWICKtype='I',reaction='Infection',time=4.403343008073078]:0.18495238495119537)[NEWICKtype='I',reaction='Infection',time=4.218390623121882]:0.8011896429483154,((77[NEWICKtype='I',reaction='Recovery',time=10.212875593577674]:3.1958770881049485,91[NEWICKtype='I',reaction='Sampling',time=7.338316663181152]:0.32131815770842564)[NEWICKtype='I',reaction='Infection',time=7.016998505472726]:1.1736030209546113,(((23[NEWICKtype='I',reaction='Recovery',time=9.725235721911726]:1.0244158012059774,35[NEWICKtype='I',reaction='Sampling',time=9.25809420796223]:0.5572742872564813)[NEWICKtype='I',reaction='Infection',time=8.700819920705749]:0.34629467868753494,((11[NEWICKtype='I',reaction='Recovery',time=10.431425785205324]:0.11660593748062809,98[NEWICKtype='I',reaction='Recovery',time=10.60324576921562]:0.2884259214909246)[NEWICKtype='I',reaction='Infection',time=10.314819847724696]:0.5460000646040939,30[NEWICKtype='I',reaction='Recovery',time=11.346861355563163]:1.5780415724425616)[NEWICKtype='I',reaction='Infection',time=9.768819783120602]:1.4142945411023877)[NEWICKtype='I',reaction='Infection',time=8.354525242018214]:1.7108614235922541,52[NEWICKtype='I',reaction='Recovery',time=8.25997757066814]:1.6163137522421795)[NEWICKtype='I',reaction='Infection',time=6.64366381842596]:0.8002683339078454)[NEWICKtype='I',reaction='Infection',time=5.843395484518115]:2.4261945043445476)[NEWICKtype='I',reaction='Infection',time=3.417200980173567]:0.06711068232398532,(55[NEWICKtype='I',reaction='Recovery',time=11.702984471804339]:4.781011530516491,17[NEWICKtype='I',reaction='Sampling',time=12.780853030716692]:5.858880089428844)[NEWICKtype='I',reaction='Infection',time=6.921972941287848]:3.5718826434382662)[NEWICKtype='I',reaction='Infection',time=3.3500902978495817]:0.46762149432250366,(96[NEWICKtype='I',reaction='Recovery',time=4.025126602568243]:0.38974046220477154,10[NEWICKtype='I',reaction='Sampling',time=4.885036824926575]:1.2496506845631044)[NEWICKtype='I',reaction='Infection',time=3.635386140363471]:0.7529173368363931)[NEWICKtype='I',reaction='Infection',time=2.882468803527078]:0.8962121252404045,(((74[NEWICKtype='I',reaction='Recovery',time=6.399619436474946]:1.8837346791718197,(((19[NEWICKtype='I',reaction='Recovery',time=5.50147661649218]:0.2571526011271228,(64[NEWICKtype='I',reaction='Recovery',time=8.275598368417768]:2.7887723990200106,59[NEWICKtype='I',reaction='Sampling',time=5.604938087831364]:0.11811211843360692)[NEWICKtype='I',reaction='Infection',time=5.486825969397757]:0.2425019540326998)[NEWICKtype='I',reaction='Infection',time=5.244324015365057]:0.15061395139649658,((62[NEWICKtype='I',reaction='Recovery',time=7.636911502265513]:0.1721368979327691,(73[NEWICKtype='I',reaction='Recovery',time=11.023702445807533]:0.7851579663529034,44[NEWICKtype='I',reaction='Recovery',time=11.808531275623597]:1.5699867961689673)[NEWICKtype='I',reaction='Infection',time=10.23854447945463]:2.7737698751218858)[NEWICKtype='I',reaction='Infection',time=7.464774604332744]:2.0226379858191086,((78[NEWICKtype='I',reaction='Recovery',time=8.815739150882727]:1.2810421528437788,(66[NEWICKtype='I',reaction='Sampling',time=8.912621113165356]:0.4800746504106126,61[NEWICKtype='I',reaction='Recovery',time=9.437281560759322]:1.0047350980045788)[NEWICKtype='I',reaction='Infection',time=8.432546462754743]:0.8978494647157946)[NEWICKtype='I',reaction='Infection',time=7.534696998038949]:1.779146624408618,((28[NEWICKtype='I',reaction='Recovery',time=9.140108896271325]:0.644896407049,46[NEWICKtype='I',reaction='Recovery',time=10.470255937537148]:1.975043448314823)[NEWICKtype='I',reaction='Infection',time=8.495212489222325]:2.25234871889561,72[NEWICKtype='I',reaction='Sampling',time=6.893199011433314]:0.6503352411065988)[NEWICKtype='I',reaction='Infection',time=6.242863770326715]:0.4873133966963845)[NEWICKtype='I',reaction='Infection',time=5.755550373630331]:0.3134137551166951)[NEWICKtype='I',reaction='Infection',time=5.4421366185136355]:0.34842655454507465)[NEWICKtype='I',reaction='Infection',time=5.093710063968561]:0.08142935346869962,(31[NEWICKtype='I',reaction='Recovery',time=6.5004079345413155]:0.8201285598215611,(49[NEWICKtype='I',reaction='Recovery',time=6.5973806990579735]:0.7732538485247655,99[NEWICKtype='I',reaction='Recovery',time=5.977427277606065]:0.15330042707285685)[NEWICKtype='I',reaction='Infection',time=5.824126850533208]:0.14384747581345358)[NEWICKtype='I',reaction='Infection',time=5.680279374719754]:0.6679986642198932)[NEWICKtype='I',reaction='Infection',time=5.012280710499861]:0.4963959531967346)[NEWICKtype='I',reaction='Infection',time=4.515884757303127]:0.7473295563109597,((((53[NEWICKtype='I',reaction='Recovery',time=8.382385135950484]:0.5483546354617461,(57[NEWICKtype='I',reaction='Recovery',time=10.650056775862737]:1.5352708787432974,(34[NEWICKtype='I',reaction='Recovery',time=9.935116595737636]:0.8017524086121206,(12[NEWICKtype='I',reaction='Recovery',time=12.00199750549129]:1.2182266602962883,(81[NEWICKtype='I',reaction='Recovery',time=12.728753835509927]:1.101438868727076,26[NEWICKtype='I',reaction='Recovery',time=12.36380535976325]:0.7364903929803983)[NEWICKtype='I',reaction='Infection',time=11.627314966782851]:0.8435441215878487)[NEWICKtype='I',reaction='Infection',time=10.783770845195003]:1.6504066580694872)[NEWICKtype='I',reaction='Infection',time=9.133364187125515]:0.018578290006075804)[NEWICKtype='I',reaction='Infection',time=9.11478589711944]:1.2807553966307017)[NEWICKtype='I',reaction='Infection',time=7.834030500488738]:1.0690773537474865,(((16[NEWICKtype='I',reaction='Sampling',time=9.710684392309021]:0.5523748816782845,(38[NEWICKtype='I',reaction='Sampling',time=9.788899411691432]:0.5363608991026343,(84[NEWICKtype='I',reaction='Recovery',time=11.591543472795488]:2.1353573546402913,71[NEWICKtype='I',reaction='Recovery',time=12.333475174143397]:2.8772890559882)[NEWICKtype='I',reaction='Infection',time=9.456186118155196]:0.20364760556639894)[NEWICKtype='I',reaction='Infection',time=9.252538512588798]:0.09422900195806072)[NEWICKtype='I',reaction='Infection',time=9.158309510630737]:0.5775522490272174,68[NEWICKtype='I',reaction='Recovery',time=9.32998830271614]:0.7492310411126208)[NEWICKtype='I',reaction='Infection',time=8.58075726160352]:0.9073952942408914,((89[NEWICKtype='I',reaction='Recovery',time=9.339375450369467]:0.3948696330409902,((13[NEWICKtype='I',reaction='Recovery',time=11.772614329761426]:2.5675759584117337,54[NEWICKtype='I',reaction='Recovery',time=12.016514433786762]:2.8114760624370696)[NEWICKtype='I',reaction='Infection',time=9.205038371349692]:0.062178166339046825,15[NEWICKtype='I',reaction='Recovery',time=12.319203551905572]:3.1763433468949263)[NEWICKtype='I',reaction='Infection',time=9.142860205010646]:0.19835438768216918)[NEWICKtype='I',reaction='Infection',time=8.944505817328476]:0.020153614641559514,2[NEWICKtype='I',reaction='Recovery',time=10.718063937967466]:1.7937117352805494)[NEWICKtype='I',reaction='Infection',time=8.924352202686917]:1.2509902353242888)[NEWICKtype='I',reaction='Infection',time=7.673361967362628]:0.9084088206213767)[NEWICKtype='I',reaction='Infection',time=6.764953146741251]:0.1739445593663662,69[NEWICKtype='I',reaction='Recovery',time=7.117028472142579]:0.5260198847676936)[NEWICKtype='I',reaction='Infection',time=6.591008587374885]:2.4114897397231223,(67[NEWICKtype='I',reaction='Recovery',time=5.055151707883041]:0.5859452551639297,40[NEWICKtype='I',reaction='Recovery',time=4.678843163551759]:0.20963671083264845)[NEWICKtype='I',reaction='Infection',time=4.469206452719111]:0.28968760506734803)[NEWICKtype='I',reaction='Infection',time=4.179518847651763]:0.41096364665959584)[NEWICKtype='I',reaction='Infection',time=3.768555200992167]:0.7095913121446329,((93[NEWICKtype='I',reaction='Recovery',time=7.89120012219306]:2.416239018964405,(86[NEWICKtype='I',reaction='Recovery',time=11.140793353371588]:0.6578700638218109,75[NEWICKtype='I',reaction='Sampling',time=11.223958565900139]:0.7410352763503614)[NEWICKtype='I',reaction='Infection',time=10.482923289549777]:5.007962186321123)[NEWICKtype='I',reaction='Infection',time=5.4749611032286545]:0.006195969720623751,((3[NEWICKtype='I',reaction='Recovery',time=8.025086191152162]:0.08311785285125595,20[NEWICKtype='I',reaction='Recovery',time=10.498544244068833]:2.556575905767927)[NEWICKtype='I',reaction='Infection',time=7.9419683383009065]:0.5270749757736688,45[NEWICKtype='I',reaction='Recovery',time=10.017020879830847]:2.602127517303609)[NEWICKtype='I',reaction='Infection',time=7.414893362527238]:1.946128229019207)[NEWICKtype='I',reaction='Infection',time=5.468765133508031]:2.4098012446604966)[NEWICKtype='I',reaction='Infection',time=3.058963888847534]:1.0727072105608606)[NEWICKtype='I',reaction='Infection',time=1.9862566782866735]:0.8728856404660741,6[NEWICKtype='I',reaction='Recovery',time=3.2379611171174685]:2.124590079296869)[NEWICKtype='I',reaction='Infection',time=1.1133710378205994]:0.4204198279390081,50[NEWICKtype='I',reaction='Recovery',time=3.3388126156559146]:2.6458614057743235)[NEWICKtype='I',reaction='Infection',time=0.6929512098815913]:0.39748021816855017)[NEWICKtype='I',reaction='Infection',time=0.29547099171304114]:0.29547099171304114;";
    public static final double TRUE_ORIGIN = 12.7808530307;
    public static final int TRUE_NS0 = 999;
    public static final double TRUE_BETA = 0.00075;
    public static final int Nt = 1001;

    public static void main(String[] args) throws Exception {

        TreeParser tree = new TreeParser(TRUE_TREE, false);

        RealParameter origin = new RealParameter(new Double[]{TRUE_ORIGIN});

        TreeIntervals intervals = new TreeIntervals(tree);

        RealParameter n_S0 = new RealParameter(new Double[] {(double)TRUE_NS0});
        RealParameter beta = new RealParameter(new Double[] {TRUE_BETA});

        // read first argument as the output file name, otherwise default to 'likelihoodCurveForStochasticSIR.txt'
        String outfileName = args.length > 0 ? args[0] : "likelihoodCurveForStochasticSIR_" + Nt + ".txt";
        System.out.println("Writing output to '" + outfileName + "'");

        // read second argument as the minimum ensemble size if present, otherwise default to 1000
        int minEnsembleSize = args.length > 1 ? Integer.parseInt(args[1]) : 1000;
        System.out.println("Using minimum ensemble size of " + minEnsembleSize);

        // read second argument as the minimum number of successful trajectories if present, otherwise default to 1000
        int minTrajSuccess = args.length > 2 ? Integer.parseInt(args[2]) : 1000;
        System.out.println("Using minimum number of successful trajectories of " + minTrajSuccess);

        // read third argument as number of ensembles size if present, otherwise default to 10
        int numberEnsemblesPerStep = args.length > 3 ? Integer.parseInt(args[3]) : 10;
        System.out.println("Number of ensembles per step " + numberEnsemblesPerStep);

        PrintWriter writer = new PrintWriter(new FileWriter(outfileName));

        writer.print("gamma\tlogP");
        System.out.print("gamma\tlogP");

        for (int i = 0; i < numberEnsemblesPerStep; i++) {
            writer.print("\tensemble_" + i);
            System.out.print("\tensemble_" + i);
        }
        for (int i = 0; i < numberEnsemblesPerStep; i++) {
            writer.print("\tensembleSize_" + i);
            System.out.print("\tensembleSize_" + i);
        }

        writer.println();
        System.out.println();

        for (double g = 0.10; g <= 0.80; g += 0.025) {

            RealParameter gamma = new RealParameter(new Double[]{g});

            writer.print(g );
            System.out.print(g);

            double[] logP = new double[numberEnsemblesPerStep];
            int[] ensembleSize = new int[numberEnsemblesPerStep];
            double maxLogP = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < numberEnsemblesPerStep; i++) {
                StochasticSIR ssir = new StochasticSIR(n_S0, beta, gamma, origin, Nt, Nt);
                StochasticSIRPopulationFunction ssirPopFun = new StochasticSIRPopulationFunction(ssir);

                StochasticCoalescent c = new StochasticCoalescent(intervals, ssirPopFun, minEnsembleSize, minTrajSuccess);

                logP[i] = c.calculateLogP();
                ensembleSize[i] = c.getLastEnsembleSize();
                if (logP[i] > maxLogP) maxLogP = logP[i];
            }

            double sumP = 0.0;
            for (int i = 0; i < numberEnsemblesPerStep; i++) {
                sumP += Math.exp(logP[i] - maxLogP);
            }
            double meanLogP = Math.log(sumP/numberEnsemblesPerStep) + maxLogP;

            writer.print("\t" + meanLogP);
            System.out.print("\t" + meanLogP);
            for (int i = 0; i < numberEnsemblesPerStep; i++) {
                writer.print("\t" + logP[i]);
                System.out.print("\t" + logP[i]);
            }

            for (int i = 0; i < numberEnsemblesPerStep; i++) {
                writer.print("\t" + ensembleSize[i]);
                System.out.print("\t" + ensembleSize[i]);
            }

            writer.println();
            System.out.println();
            writer.flush();
        }
        writer.close();
    }
}
