from covigator.database.model import VariantType

VARIANT_TYPE_COLOR_MAP = {
    VariantType.SNV.name: "#8da0cb",
    VariantType.INSERTION.name: "#fc8d62",
    VariantType.DELETION.name: "#66c2a5",
}