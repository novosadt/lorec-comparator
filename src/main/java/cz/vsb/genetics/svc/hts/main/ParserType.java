package cz.vsb.genetics.svc.hts.main;

import java.util.HashMap;
import java.util.Map;

public enum ParserType {
    BIONANO("b"),
    ANNOT_SV("a"),
    SAMPLOT("s"),
    VCF_LONGRANGER("vl"),
    VCF_SNIFFLES("vs"),
    VCF_MANTA("vm"),
    VCF_ICLR("vi")
    
    ;

    public final String value;

    private static final Map<String, ParserType> map = new HashMap<>();

    static
    {
        for (ParserType item : ParserType.values())
            map.put(item.value, item);
    }

    ParserType(String value) {
        this.value = value;
    }

    public static ParserType of(String value) {
        return map.get(value);
    }
}
