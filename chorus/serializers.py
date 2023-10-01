import os

from rest_flex_fields import FlexFieldsModelSerializer
from rest_framework import serializers

from chorus.models import Protein, Variant, ChorusSession


class ProteinSerializer(FlexFieldsModelSerializer):
    variants = serializers.SerializerMethodField()

    def get_variants(self, obj):
        data = []
        for i in obj.variant_set.all():
            data.append({
                "position": i.position,
                "original": i.original,
                "mutated": i.mutated,
                "score": i.score,
                "pathogenicity": i.pathogenicity
            })
        return data

    class Meta:
        model = Protein
        fields = ["id", "name", "description", "variants"]


class VariantSerializer(FlexFieldsModelSerializer):
    protein = serializers.SerializerMethodField()

    def get_protein(self, obj):
        return obj.protein.name

    class Meta:
        model = Variant
        fields = ["id", "protein", "position", "original", "mutated", "score", "pathogenicity"]


class ChorusSessionSerializer(FlexFieldsModelSerializer):
    file = serializers.SerializerMethodField()

    def get_file(self, record):
        _, filename = os.path.split(record.file.name)
        return filename

    class Meta:
        model = ChorusSession
        fields = ["id", "user", "link_id", "description", "file", "created_at", "updated_at"]
        lookup_field = "link_id"
