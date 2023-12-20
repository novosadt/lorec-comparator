/*
 * om-hts-svc: Optical Mapping and High-Throughput-Sequencing Variant Comparator
 *
 * Application for comparison of structural variants found by optical mapping technology (Bionano Genomics)
 * with AnnotSV and Samplot analysis of 3rd generation sequencing technologies 10xGenomics, Oxford Nanopore Technologies and Pacbio.
 *
 *
 * Copyright (C) 2022  Tomas Novosad
 * VSB-TUO, Faculty of Electrical Engineering and Computer Science
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package cz.vsb.genetics.svc.hts.main;

import cz.vsb.genetics.common.ChromosomeRegion;
import cz.vsb.genetics.ngs.sv.AnnotSvTsvParser;
import cz.vsb.genetics.ngs.sv.GenericSvVcfParser;
import cz.vsb.genetics.ngs.sv.SamplotCsvParser;
import cz.vsb.genetics.om.sv.BionanoPipelineResultParser;
import cz.vsb.genetics.sv.MultipleSvComparator;
import cz.vsb.genetics.sv.StructuralVariantType;
import cz.vsb.genetics.sv.SvResultParser;
import org.apache.commons.cli.*;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

public class BionanoHtsSvComparator {
    private static final Logger log = LoggerFactory.getLogger(BionanoHtsSvComparator.class);

    private static final String ARG_BIONANO_INPUT = "bionano_input";
    private static final String ARG_ANNOTSV_INPUT = "annotsv_input";
    private static final String ARG_SAMPLOT_INPUT = "samplot_input";
    private static final String ARG_VCF_LONGRANGER_INPUT = "vcf_longranger_input";
    private static final String ARG_VCF_SNIFFLES_INPUT = "vcf_sniffles_input";
    private static final String ARG_VCF_MANTA_INPUT = "vcf_manta_input";
    private static final String ARG_VCF_ICLR_INPUT = "vcf_iclr_input";
    private static final String ARG_MAIN_INPUT = "main_input";
    private static final String ARG_VCF_FILTER_PASS = "vcf_filter_pass";
    private static final String ARG_VARIANT_TYPE = "variant_type";
    private static final String ARG_DISTANCE_VARIANCE = "distance_variance";
    private static final String ARG_INTERSECTION_VARIANCE = "intersection_variance";
    private static final String ARG_MINIMAL_PROPORTION = "minimal_proportion";
    private static final String ARG_GENE_INTERSECTION = "gene_intersection";
    private static final String ARG_PREFER_BASE_SVTYPE = "prefer_base_svtype";
    private static final String ARG_STATISTICS_OUTPUT = "statistics_output";
    private static final String ARG_DISTANCE_VARIANCE_STATISTICS = "distance_variance_statistics";
    private static final String ARG_INTERSECTION_VARIANCE_STATISTICS = "intersection_variance_statistics";
    private static final String ARG_REGION_FILTER_FILE = "region_filter_file";
    private static final String ARG_OUTPUT = "output";

    private SvResultParser mainParser;
    private final List<SvResultParser> otherParsers = new ArrayList<>();

    public static void main(String[] args) {
        try {
            BionanoHtsSvComparator comparator = new BionanoHtsSvComparator();
            comparator.compareVariants(args);
        }
        catch (Exception e) {
            System.err.println("Error occurred: " + e.getMessage());
            e.printStackTrace();
        }
    }

    public void compareVariants(String[] args) throws Exception {
        CommandLine cmd = getCommandLine(args);

        Integer distanceVariance = cmd.hasOption(ARG_DISTANCE_VARIANCE) ? Integer.valueOf(cmd.getOptionValue(ARG_DISTANCE_VARIANCE)) : null;
        Double intersectionVariance = cmd.hasOption(ARG_INTERSECTION_VARIANCE) ? new Double(cmd.getOptionValue(ARG_INTERSECTION_VARIANCE)) : null;
        Double minimalProportion = cmd.hasOption(ARG_MINIMAL_PROPORTION) ? new Double(cmd.getOptionValue(ARG_MINIMAL_PROPORTION)) : null;
        Set<StructuralVariantType> variantType = cmd.hasOption(ARG_VARIANT_TYPE) ? StructuralVariantType.getSvTypes(cmd.getOptionValue(ARG_VARIANT_TYPE)) : null;
        boolean onlyCommonGeneVariants = cmd.hasOption(ARG_GENE_INTERSECTION);
        String statsOutput = cmd.getOptionValue(ARG_STATISTICS_OUTPUT);
        String[] distanceVarianceStatsCounts = cmd.hasOption(ARG_DISTANCE_VARIANCE_STATISTICS) ? cmd.getOptionValue(ARG_DISTANCE_VARIANCE_STATISTICS).split(";") : null;
        String[] intersectionVarianceStatsThreshold = cmd.hasOption(ARG_INTERSECTION_VARIANCE_STATISTICS) ? cmd.getOptionValue(ARG_INTERSECTION_VARIANCE_STATISTICS).split(";") : null;

        boolean calculateStats = statsOutput != null &&
                ((distanceVarianceStatsCounts != null && distanceVarianceStatsCounts.length > 0) ||
                 (intersectionVarianceStatsThreshold != null && intersectionVarianceStatsThreshold.length > 0));

        initParsers(cmd);

        MultipleSvComparator svComparator = new MultipleSvComparator();
        svComparator.setOnlyCommonGenes(onlyCommonGeneVariants);
        svComparator.setDistanceVarianceThreshold(distanceVariance);
        svComparator.setIntersectionVarianceThreshold(intersectionVariance);
        svComparator.setVariantTypes(variantType);
        svComparator.setMinimalProportion(minimalProportion);
        svComparator.setExcludedRegions(getExcludedRegions(cmd.getOptionValue(ARG_REGION_FILTER_FILE)));

        if (calculateStats) {
            svComparator.setCalculateStructuralVariantStats(true);

            if (distanceVarianceStatsCounts != null && distanceVarianceStatsCounts.length > 0)
                svComparator.setDistanceVarianceBasesCounts(Arrays.stream(distanceVarianceStatsCounts).mapToInt(Integer::parseInt).toArray());

            if (intersectionVarianceStatsThreshold != null && intersectionVarianceStatsThreshold.length > 0)
                svComparator.setIntersectionVarianceThresholds(Arrays.stream(intersectionVarianceStatsThreshold).mapToDouble(Double::parseDouble).toArray());
        }

        svComparator.compareStructuralVariants(mainParser, otherParsers, cmd.getOptionValue(ARG_OUTPUT));

        printStructuralVariants(mainParser, otherParsers);

        if (calculateStats)
            svComparator.saveStructuralVariantStats(statsOutput);
    }

    private CommandLine getCommandLine(String[] args) {
        Options options = new Options();

        Option bionanoInput = new Option("b", ARG_BIONANO_INPUT, true, "bionano pipeline result file path (smap)");
        bionanoInput.setRequired(true);
        bionanoInput.setArgName("smap file");
        bionanoInput.setType(String.class);
        options.addOption(bionanoInput);

        Option annotsvInput = new Option("a", ARG_ANNOTSV_INPUT, true, "annotsv tsv file paths delimited by semicolon");
        annotsvInput.setArgName("tsv file");
        annotsvInput.setType(String.class);
        options.addOption(annotsvInput);

        Option samplotVariants = new Option("s", ARG_SAMPLOT_INPUT, true, "samplot csv variants file paths delimited by semicolon");
        samplotVariants.setArgName("csv file");
        samplotVariants.setType(String.class);
        options.addOption(samplotVariants);

        Option vcfLongrangerInput = new Option("vl", ARG_VCF_LONGRANGER_INPUT, true, "longranger vcf variants file paths delimited by semicolon");
        vcfLongrangerInput.setArgName("vcf file");
        vcfLongrangerInput.setType(String.class);
        options.addOption(vcfLongrangerInput);

        Option vcfSnifflesInput = new Option("vs", ARG_VCF_SNIFFLES_INPUT, true, "sniffles vcf variants file paths delimited by semicolon");
        vcfSnifflesInput.setArgName("vcf file");
        vcfSnifflesInput.setType(String.class);
        options.addOption(vcfSnifflesInput);

        Option vcfMantaInput = new Option("vm", ARG_VCF_MANTA_INPUT, true, "manta vcf variants file paths delimited by semicolon");
        vcfMantaInput.setArgName("vcf file");
        vcfMantaInput.setType(String.class);
        options.addOption(vcfMantaInput);

        Option vcfIclrInput = new Option("vi", ARG_VCF_ICLR_INPUT, true, "Illumina Dragen ICLR wgs vcf variants file paths delimited by semicolon");
        vcfIclrInput.setArgName("vcf file");
        vcfIclrInput.setType(String.class);
        options.addOption(vcfIclrInput);

        Option mainInput = new Option("mi", ARG_MAIN_INPUT, true, "Main variant file path used to determine main technology and input between other inputs.");
        mainInput.setArgName("file path");
        mainInput.setType(String.class);
        options.addOption(mainInput);

        Option vcfFilterPass = new Option("vfp", ARG_VCF_FILTER_PASS, false, "Process only structural variants with filter value PASS");
        options.addOption(vcfFilterPass);

        Option geneIntersection = new Option("g", ARG_GENE_INTERSECTION, false, "select only variants with common genes (default false)");
        options.addOption(geneIntersection);

        Option svType = new Option("svt", ARG_PREFER_BASE_SVTYPE, false, "whether to prefer base variant type (SVTYPE) in case of BND and 10x/TELL-Seq (default false i.e. SVTYPE2)");
        options.addOption(svType);

        Option variantType = new Option("t", ARG_VARIANT_TYPE, true, "variant type filter, any combination of [BND,CNV,DEL,INS,DUP,INV,UNK], delimited by semicolon");
        variantType.setType(String.class);
        variantType.setArgName("sv types");
        options.addOption(variantType);

        Option distanceVariance = new Option("d", ARG_DISTANCE_VARIANCE, true, "distance variance filter - number of bases difference between variant from NGS and OM");
        distanceVariance.setType(Integer.class);
        distanceVariance.setArgName("number of bases");
        distanceVariance.setRequired(false);
        options.addOption(distanceVariance);

        Option intersectionVariance = new Option("i", ARG_INTERSECTION_VARIANCE, true, "intersection variance filter - threshold difference between variant from NGS and OM");
        intersectionVariance.setType(Integer.class);
        intersectionVariance.setArgName("threshold");
        intersectionVariance.setRequired(false);
        options.addOption(intersectionVariance);

        Option minimalProportion = new Option("mp", ARG_MINIMAL_PROPORTION, true, "minimal proportion filter - minimal proportion of target variant within query variant (0.0 - 1.0)");
        minimalProportion.setType(Double.class);
        minimalProportion.setArgName("number");
        minimalProportion.setRequired(false);
        options.addOption(minimalProportion);

        Option distanceVarianceStats = new Option("dvs", ARG_DISTANCE_VARIANCE_STATISTICS, true, "distance variance statistics - bases counts delimited by semicolon (e.g. 10000;50000;100000)");
        distanceVarianceStats.setType(Integer.class);
        distanceVarianceStats.setArgName("bases counts");
        distanceVarianceStats.setRequired(false);
        options.addOption(distanceVarianceStats);

        Option intersectionVarianceStats = new Option("ivs", ARG_INTERSECTION_VARIANCE_STATISTICS, true, "intersection variance statistics - thresholds delimited by semicolon (e.g. 0.1;0.3;0.5)");
        intersectionVarianceStats.setType(Integer.class);
        intersectionVarianceStats.setArgName("thresholds");
        intersectionVarianceStats.setRequired(false);
        options.addOption(intersectionVarianceStats);

        Option statsOutput = new Option("so", ARG_STATISTICS_OUTPUT, true, "output structural variants statistics file");
        statsOutput.setArgName("csv file");
        statsOutput.setType(String.class);
        options.addOption(statsOutput);

        Option regionFilterFile = new Option("rff", ARG_REGION_FILTER_FILE, true, "List of regions to be excluded from analysis (bed format, tab separated)");
        regionFilterFile.setArgName("bed file");
        regionFilterFile.setType(String.class);
        options.addOption(regionFilterFile);

        Option output = new Option("o", ARG_OUTPUT, true, "output result file");
        output.setRequired(true);
        output.setArgName("csv file");
        output.setType(String.class);
        options.addOption(output);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println("\nSVC - Bionano Genomics (OM) and High Throughput Sequencing (HTS) Structural Variant Comparator, v" + BionanoHtsSvComparator.version() + "\n");
            System.out.println(e.getMessage());
            System.out.println();
            formatter.printHelp(
                    300,
                    "\njava -jar om-samplot-svc.jar ",
                    "\noptions:",
                    options,
                    "\nTomas Novosad, VSB-TU Ostrava, 2023" +
                            "\nFEI, Department of Computer Science" +
                            "\nVersion: " + version() +
                            "\nLicense: GPL-3.0-only ");

            System.exit(1);
        }

        return cmd;
    }

    private void initParsers(CommandLine cmd) throws Exception {
        boolean preferBaseSvType = cmd.hasOption(ARG_PREFER_BASE_SVTYPE);
        boolean vcfFilterPass = cmd.hasOption(ARG_VCF_FILTER_PASS);
        String mainInput = cmd.hasOption(ARG_MAIN_INPUT) ? cmd.getOptionValue(ARG_MAIN_INPUT) : cmd.getOptionValue(ARG_BIONANO_INPUT);

        if (StringUtils.isBlank(mainInput) && !cmd.hasOption(ARG_BIONANO_INPUT)) {
            System.out.println("Cannon determine main input source. Either -mi parameter or Bionano smap input must be specified.");
            System.exit(1);
        }

        if (cmd.hasOption(ARG_BIONANO_INPUT)) {
            String input = cmd.getOptionValue(ARG_BIONANO_INPUT);
            SvResultParser bionanoParser = new BionanoPipelineResultParser("bionano");
            bionanoParser.setRemoveDuplicateVariants(true);
            bionanoParser.parseResultFile(input, "[,\t]");
            addParser(bionanoParser, input.equals(mainInput));
        }

        if (cmd.hasOption(ARG_ANNOTSV_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_ANNOTSV_INPUT).split(";");

            for (String input : inputs) {
                SvResultParser annotsvParser = new AnnotSvTsvParser("annotsv_" + getParserNameSuffix(input), preferBaseSvType);
                annotsvParser.setRemoveDuplicateVariants(true);
                annotsvParser.parseResultFile(input, "\t");
                addParser(annotsvParser, input.equals(mainInput));
            }
        }

        if (cmd.hasOption(ARG_SAMPLOT_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_SAMPLOT_INPUT).split(";");

            for (String input : inputs) {
                SvResultParser samplotParser = new SamplotCsvParser("samplot_" + getParserNameSuffix(input));
                samplotParser.setRemoveDuplicateVariants(true);
                samplotParser.parseResultFile(input, "\t");
                addParser(samplotParser, input.equals(mainInput));
            }
        }

        if (cmd.hasOption(ARG_VCF_LONGRANGER_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_VCF_LONGRANGER_INPUT).split(";");
            addVcfOtherParser("vcf-longranger_", inputs, vcfFilterPass, preferBaseSvType, mainInput);
        }

        if (cmd.hasOption(ARG_VCF_SNIFFLES_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_VCF_SNIFFLES_INPUT).split(";");
            addVcfOtherParser("vcf-sniffles_", inputs, vcfFilterPass, preferBaseSvType, mainInput);
        }

        if (cmd.hasOption(ARG_VCF_MANTA_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_VCF_MANTA_INPUT).split(";");
            addVcfOtherParser("vcf-manta_", inputs, vcfFilterPass, preferBaseSvType, mainInput);
        }

        if (cmd.hasOption(ARG_VCF_ICLR_INPUT)) {
            String[] inputs = cmd.getOptionValue(ARG_VCF_ICLR_INPUT).split(";");
            addVcfOtherParser("vcf-iclr_", inputs, vcfFilterPass, preferBaseSvType, mainInput);
        }
    }

    private void addParser(SvResultParser parser, boolean main) {
        if (main)
            mainParser = parser;
        else
            otherParsers.add(parser);
    }

    private void addVcfOtherParser(String namePrefix, String[] inputs, boolean vcfFilterPass, boolean preferBaseSvType, String mainInput) throws Exception {
        for (String input : inputs) {
            GenericSvVcfParser vcfParser = new GenericSvVcfParser(namePrefix + getParserNameSuffix(input));
            vcfParser.setOnlyFilterPass(vcfFilterPass);
            vcfParser.setPreferBaseSvType(preferBaseSvType);
            vcfParser.setRemoveDuplicateVariants(true);
            vcfParser.parseResultFile(input, "\t");
            addParser(vcfParser, input.equals(mainInput));
        }

    }

    private String getParserNameSuffix(String name) {
        return FilenameUtils.removeExtension(new File(name).getName());
    }

    private void printStructuralVariants(SvResultParser mainParser, List<SvResultParser> otherParsers) {
        mainParser.printStructuralVariantStats();

        for (SvResultParser otherParser : otherParsers)
            otherParser.printStructuralVariantStats();
    }

    private static String version() {
        final Properties properties = new Properties();

        try {
            properties.load(BionanoHtsSvComparator.class.getClassLoader().getResourceAsStream("project.properties"));
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        return properties.getProperty("version");
    }

    private List<ChromosomeRegion> getExcludedRegions(String regionFilterFile) throws Exception {
        if (StringUtils.isBlank(regionFilterFile))
            return null;

        List<ChromosomeRegion> regions = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(regionFilterFile))) {
            String line;
            while((line = reader.readLine()) != null) {
                ChromosomeRegion region = ChromosomeRegion.valueOf(line, "\t");

                if (region != null)
                    regions.add(region);
            }
        }
        return regions;
    }
}
