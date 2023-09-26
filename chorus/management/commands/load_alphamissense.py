import re

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from chorus.models import Protein, Variant

class Command(BaseCommand):
    """
    A command to parse alphamissense tabulated txt file and store it in database using django model Protein for uniprot id and variant for the actual missense position, mutation and score
    """
    help = 'Parse alphamissense tabulated txt file and store it in database using django model Protein for uniprot id and variant for the actual missense position, mutation and score'

    def add_arguments(self, parser):
        parser.add_argument('file_path', type=str, help='Path to the alphamissense tabulated txt file to be processed')

    def handle(self, *args, **options):
        file_path = options['file_path']
        pattern = re.compile(r"([A-Z]+)(\d+)([A-Z]+)")
        with open(file_path, 'r') as f:
            temp_data = []
            current_protein = None
            for _ in range(4):
                next(f)
            for line in f:
                line = line.strip().split('\t')
                protein_name = line[0]
                search_result = pattern.search(line[1])
                if search_result:
                    score = line[2]
                    position = int(search_result.group(2))
                    original = search_result.group(1)
                    mutated = search_result.group(3)
                    pathogenicity = line[3]
                    if current_protein is None:
                        protein = Protein.objects.get_or_create(name=protein_name, description="")
                        if type(protein) == tuple:
                            current_protein = protein[0]
                        else:
                            current_protein = protein
                    if current_protein.name != protein_name:
                        with transaction.atomic():
                            for data in temp_data:
                                Variant.objects.create(protein=current_protein, position=data[1], original=data[2], mutated=data[3], score=data[4], pathogenicity=data[5])
                            temp_data = []
                            protein = Protein.objects.get_or_create(name=protein_name, description="")
                            if type(protein) == tuple:
                                current_protein = protein[0]
                            else:
                                current_protein = protein
                    temp_data.append((protein_name, position, original, mutated, score, pathogenicity))
            with transaction.atomic():
                for data in temp_data:
                    Variant.objects.create(protein=current_protein, position=data[1], original=data[2], mutated=data[3], score=data[4], pathogenicity=data[5])
        self.stdout.write(self.style.SUCCESS('Successfully parsed and stored alphamissense data in database'))