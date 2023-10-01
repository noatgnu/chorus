import io

import pandas as pd
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from filters.mixins import FiltersMixin
from rest_framework import permissions, filters, status
from rest_framework.decorators import action
from rest_framework.parsers import MultiPartParser, JSONParser
from rest_framework.response import Response
from rest_framework.viewsets import ModelViewSet

from chorus.filter_schema import protein_query_schema, chorus_session_query_schema
from chorus.models import Protein, Variant, ChorusSession
from chorus.serializers import ProteinSerializer, VariantSerializer, ChorusSessionSerializer

from uniprotparser.betaparser import UniprotParser

from chorus.utils import extract_functional_domain
from django.core.files.base import File as djangoFile

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
        result = list(parser.parse([protein.name]))
        if len(result) > 0:
            df = pd.read_csv(io.StringIO(result[0]), sep="\t")
            df = df.apply(lambda x: extract_functional_domain(x, "Domain [FT]"), axis=1)
            df.fillna("", inplace=True)
            return Response(df.to_dict(orient="records")[0])
        else:
            return Response(status=404)


class VariantViewSets(FiltersMixin, ModelViewSet):
    queryset = Variant.objects.all()
    serializer_class = VariantSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    filter_backends = [filters.OrderingFilter,]
    ordering_fields = ["protein", "position", "original", "mutated", "pathogenicity"]
    filter_mappings = {
        "protein": "protein__name__exact",
        "position": "position__exact",
        "original": "original__exact",
        "mutated": "mutated__exact",
        "pathogenicity": "pathogenicity__exact",
    }
    filter_validation_schema = protein_query_schema

    def get_queryset(self):
        return self.queryset


class ChorusSessionViewSets(ModelViewSet):
    queryset = ChorusSession.objects.all()
    parser_classes = [MultiPartParser, JSONParser]
    serializer_class = ChorusSessionSerializer
    permission_classes = (permissions.AllowAny,)
    filter_backends = [filters.OrderingFilter,]
    ordering_fields = ["created_at", "updated_at"]
    filter_mappings = {
        "user": "user__id__exact",
        "link_id": "link_id__exact",
        "id": "id__exact",
    }
    filter_validation_schema = chorus_session_query_schema

    def get_queryset(self):
        return self.queryset

    def create(self, request, *args, **kwargs):
        c = ChorusSession()
        if request.user:
            if request.user.is_authenticated:
                c.user = request.user
        c.file.save(f"{str(c.link_id)}.json",djangoFile(self.request.data['file']))
        c.description = request.data["description"]
        c.save()
        json_data = ChorusSessionSerializer(c, many=False, context={"request": request})
        return Response(json_data.data, status=201)

    @action(methods=["get"], url_path="download/", detail=True, permission_classes=[
        permissions.AllowAny
    ])
    @method_decorator(cache_page(0))
    def download(self, request, pk=None, link_id=None, token=None):
        c = self.get_object()
        return Response(data={"url": c.file.url}, status=status.HTTP_200_OK)