from rest_flex_fields import FlexFieldsModelSerializer

from chorus.models import Protein


class ProteinSerializer(FlexFieldsModelSerializer):
    class Meta:
        model = Protein
        fields = ["id", "name", "description"]


class VariantSerializer(FlexFieldsModelSerializer):
    class Meta:
        model = Protein
        fields = ["id", "protein", "position", "original", "mutated", "score", "pathogenicity"]