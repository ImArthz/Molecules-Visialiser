# Molecule Visualizer

## Introdução

O Molecule Visualizer é um aplicativo de visualização de moléculas que oferece uma interface gráfica para exibir moléculas em 2D e 3D. Utiliza a biblioteca RDKit para manipulação de moléculas e py3Dmol para visualizações 3D. Este README fornece uma visão geral das funcionalidades e da estrutura do código.

## Bibliotecas Necessárias

# Para executar o Molecule Visualizer, você precisará das seguintes bibliotecas:

- **RDKit**: Biblioteca para manipulação e visualização de moléculas químicas.
- **py3Dmol**: Biblioteca para visualização 3D de moléculas.
- **Tkinter**: Biblioteca padrão do Python para criar interfaces gráficas.
- **Pillow**: Biblioteca para manipulação de imagens.

# Para instalar essas bibliotecas, use os seguintes comandos:

```bash
pip install rdkit
pip install py3Dmol
pip install pillow
```
* Nota: O RDKit pode ser difícil de instalar diretamente via pip, pois pode precisar de compilação. Em alguns casos, é melhor usar uma distribuição como o Anaconda, que já inclui o RDKit.

# Passos para Executar o Código
* Preparação do Ambiente
* Crie um diretório de trabalho:
```bash
mkdir molecule_visualizer
cd molecule_visualizer
```
* Execute o Código
* Abra um terminal e navegue até o diretório onde você salvou o arquivo Python.
* Execute o script Python com o comando:
```bash
python main.py
```
* A GUI (Interface Gráfica do Usuário) será exibida. Nela, você poderá buscar moléculas pelo nome e visualizá-las em 2D e 3D.

##  Estrutura do Código
# 1. Importação de Módulos
```bash

import json
import os
import logging
import webbrowser
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
```
# 2. Função ensure_output_directory()
* Garante que o diretório outputs exista.

# 3. Função load_molecules()
* Carrega as informações das moléculas de um arquivo JSON.

# 4. Função visualize_molecule()
* Cria a estrutura 2D da molécula a partir da string SMILES.

# 5. Função draw_molecule()
* Desenha a molécula em 2D e salva como um arquivo PNG.

# 6. Função visualize_with_py3dmol()
* Gera uma visualização 3D da molécula e salva como um arquivo HTML.

# 7. Função cleanup_files()
* Remove arquivos antigos do diretório outputs.

# 8. Função update_molecule_list()
* Atualiza a lista de moléculas na GUI com base no arquivo JSON.

# 9. Função on_search_change()
* Filtra a lista de moléculas com base no termo de busca inserido pelo usuário.

# 10. Função visualize_molecule_gui()
* Lida com a visualização de uma molécula selecionada na GUI. Gera imagens 2D e 3D da molécula e abre a visualização 3D no navegador.

# 11. Função open_git_link()
* Abre um link para o repositório Git dos criadores.

# 12. Função open_chemistry_resource_link()
* Abre um link para uma página com uma explicação detalhada da molécula pesquisada.

# 13. Função exit_gui()
* Limpa arquivos temporários e fecha a aplicação.

# 14. Função main_gui()
* Cria e configura a interface gráfica principal, incluindo os novos botões para visualizar moléculas, sair da aplicação, acessar o repositório Git e abrir o recurso de química.

# Novidades
* Adição do Botão "Chemistry Resource"
* O aplicativo agora inclui um botão "Chemistry Resource" que abre um site renomado para uma explicação detalhada da molécula pesquisada.

## Como usar:
* Digite o nome da molécula no campo de busca.
* Selecione a molécula da lista exibida.
* Clique no botão "Chemistry Resource" para ser redirecionado a uma página com informações detalhadas sobre a molécula.
* Nota: O link para a página de química detalhada é um exemplo e pode precisar ser ajustado conforme a disponibilidade dos recursos.

## Conclusão
O Molecule Visualizer proporciona uma maneira poderosa e interativa de visualizar moléculas em 2D e 3D. Cada função foi projetada para modularizar e otimizar o processo de visualização e interação com as moléculas. Com uma interface gráfica amigável e recursos avançados de visualização, o aplicativo é uma ferramenta valiosa para pesquisadores e entusiastas da química.

## Colaboradores
* Arthur De Oliveira Mendonça

## Redes Sociais:

<div align="center"> <h3>Connect with me:</h3> <p><strong>Arthur Mendonça</strong></p> <p> <a href="https://github.com/ImArthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/github.svg" alt="ImArthz" height="30" width="40" /> </a> <a href="https://twitter.com/Im_Arthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/twitter.svg" alt="Im_Arthz" height="30" width="40" /> </a> <a href="https://api.whatsapp.com/send?phone=37988528423" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/whatsapp.svg" alt="WhatsApp" height="30" width="40" /> </a> <a href="https://discordapp.com/users/imarthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/discord.svg" alt="imarthz" height="30" width="40" /> </a> </p> </div>
