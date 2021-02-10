package abstract_url_parser;

import java.util.List;
import java.util.stream.Collectors;

import com.opencsv.ICSVParser;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVReader;
import com.linkedin.urls.detection.UrlDetector;
import com.linkedin.urls.detection.UrlDetectorOptions;
import com.opencsv.CSVParserBuilder;


public class App 
{
    public static void main( final String[] args) throws IOException {
        final CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(args[0]));
        final ICSVParser parser = new CSVParserBuilder().withSeparator('\t').build();
        final CSVReader csvReader = builder.withCSVParser(parser).build();
        final List<String[]> rows = csvReader.readAll();
        // remove header
        rows.remove(0);

        final File outf = new File(args[1]);
        if (!outf.exists()) {
            outf.createNewFile();
        }
        final FileWriter writer = new FileWriter(outf);
        writer.write("PMID\tURL_BRACKET\tURL_BRACKET_NORM\tURL_XML\tURL_XML_NORM\n");

        for (String[] entry : rows) {
            writer.write(entry[0]);
            writer.write('\t');
            // remove trailing dot
            if (entry[1].endsWith(".")){
                entry[1] = entry[1].substring(0, entry[1].length()-1);
            }
            final UrlDetector urlparser = new UrlDetector(entry[1], UrlDetectorOptions.BRACKET_MATCH);
            final List<String> orig_urls = urlparser.detect().stream().map(u -> u.getOriginalUrl())
                    .collect(Collectors.toList());
            final List<String> norm_urls = urlparser.detect().stream().map(u -> u.getFullUrlWithoutFragment())
                    .collect(Collectors.toList());
            final List<String> norm_urls2 = urlparser.detect().stream().map(u -> u.normalize().getFullUrlWithoutFragment())
                    .collect(Collectors.toList());
            writer.write(String.join("; ", orig_urls));
            writer.write("\t");
            writer.write(String.join("; ", norm_urls2));
            writer.write("\t");
            UrlDetector urlparser_xml = new UrlDetector(entry[1], UrlDetectorOptions.XML);
            List<String> orig_urls_xml = urlparser_xml.detect().stream().map(u -> u.getOriginalUrl()).collect(Collectors.toList());
            final List<String> norm_urls_xml = urlparser_xml.detect().stream().map(u -> u.getFullUrlWithoutFragment()).collect(Collectors.toList());
            final List<String> norm_urls_xml2 = urlparser_xml.detect().stream().map(u -> u.normalize().getFullUrlWithoutFragment()).collect(Collectors.toList());
            writer.write(String.join("; ", orig_urls_xml));
            writer.write("\t");
            writer.write(String.join("; ", norm_urls_xml2));
            writer.write("\n");
        }
        writer.close();
    }
}
