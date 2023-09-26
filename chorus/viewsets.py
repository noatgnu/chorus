import io

import pandas as pd
from filters.mixins import FiltersMixin
from rest_framework import permissions, filters
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.viewsets import ModelViewSet

from chorus.filter_schema import protein_query_schema
from chorus.models import Protein
from chorus.serializers import ProteinSerializer

from uniprotparser.betaparser import UniprotParser
class ProteinViewSets(FiltersMixin, ModelViewSet):
    queryset = Protein.objects.all()
    serializer_class = ProteinSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    filter_backends = [filters.OrderingFilter,]
    ordering_fields = ["name"]
    filter_mappings = {
        "name": "name__icontains",
    }
    filter_validation_schema = protein_query_schema

    def get_queryset(self):
        return self.queryset

    @action(detail=True, methods=["get"])
    def get_uniprot(self, request, pk=None):
        protein = self.get_object()
        parser = UniprotParser(columns="accession,id,ft_domain,sequence,gene_names")
        result = list(parser.parse(protein.name))
        if len(result) > 0:
            df = pd.read_csv(io.StringIO(result[0]), sep="\t")
            return Response(df.to_dict(orient="records"))
        else:
            return Response(status=404)


class VariantViewSets(FiltersMixin, ModelViewSet):
    queryset = Protein.objects.all()
    serializer_class = ProteinSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    filter_backends = [filters.OrderingFilter,]
    ordering_fields = ["protein", "position", "original", "mutated", "pathogenicity"]
    filter_mappings = {
        "name": "name__icontains",
        "protein": "protein__name__icontains",
        "position": "position__exact",
        "original": "original__exact",
        "mutated": "mutated__exact",
        "pathogenicity": "pathogenicity__exact",
    }
    filter_validation_schema = protein_query_schema

    def get_queryset(self):
        return self.queryset


