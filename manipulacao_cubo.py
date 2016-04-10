from random import randint
from os import system
from sys import argv
import pyperclip
import struct
import random

def gerar_atributos(cardinalidade, tuplas):
	col = list()

	if tuplas != cardinalidade:
		for i in range(1, (tuplas-cardinalidade)+1):
			col += [randint(1, cardinalidade)]

		for i in range(1, cardinalidade + 1):
			n = randint(0, tuplas - cardinalidade)
			col.insert(n, i)
	else:
		for i in range(1, cardinalidade + 1):
			index = randint(0, len(col))
			print('Insert value ' + str(i) + ' on index ' + str(index))
			col.insert(index, i)

	return col


def gerar_base_binaria(nome, cardinalidade, dimensoes, tuplas, tid):
	nome += '.bin'
	primeira_tupla = ''

	atributos = list()

	arquivo = open(nome, 'wb')
	arquivo.write(struct.pack('i', int(tuplas)))

	for i in range(1, dimensoes + 1):
		atributos += [gerar_atributos(cardinalidade, tuplas)]		
		arquivo.write(struct.pack('i', int(cardinalidade)))

	for t in range(0, tuplas):
		for d in range(0, dimensoes):
			atributo = atributos[d][t]
			arquivo.write(struct.pack('i', atributo))

			if t + 1 == tid:
				primeira_tupla += str(atributo) + ' '

	primeira_tupla = primeira_tupla[:-1]
	arquivo.close()

	return [nome, primeira_tupla]


def gerar_base_texto(nome, cardinalidade, dimensoes, tuplas, tid):
	nome += '.txt'
	primeira_tupla = ''

	atributos = list()

	arquivo = open(nome, 'w')

	arquivo.write(str(tuplas) + ' ')

	for i in range(0, dimensoes):
		atributos += [gerar_atributos(cardinalidade, tuplas)]
		arquivo.write(str(cardinalidade) + ' ')

	arquivo.write('\n')

	for t in range(0, tuplas):
		for d in range(0, dimensoes):
			atributo = atributos[d][t]
			arquivo.write(str(atributo) + ' ')

			if t + 1 == tid:
				primeira_tupla += str(atributo) + ' '

		arquivo.write('\n')	

	primeira_tupla = primeira_tupla[:-1]
	arquivo.close()

	return [nome, primeira_tupla]


def gerar_cubo():
	cardinalidade = argv[2]
	dimensoes = argv[3]
	tuplas = argv[4]
	tid = argv[5]
	tipo = argv[6]

	nome = 'novocubo_C' + cardinalidade + '_D' + dimensoes + '_T' + tuplas

	if tipo == 'binario':
		retorno = gerar_base_binaria(nome, int(cardinalidade), int(dimensoes), int(tuplas), int(tid))
	elif tipo == 'texto':
		retorno = gerar_base_texto(nome, int(cardinalidade), int(dimensoes), int(tuplas), int(tid))		
	else:
		print('Erro de parametrização! Abortando...')
		sys.exit(0)
	
	nome = retorno[0]
	consulta = retorno[1]

	pyperclip.copy(retorno[1])

	print('\n\nConsulta \"' + consulta + '\" copiada para a clipboard!')
	print('\n\nArquivo ' + nome + ' gerado com sucesso!')


def converter_novocubo_binario_pra_texto(p_arq_binario, p_arq_texto, dimensoes):
	arq_binario = open(p_arq_binario, 'rb')
	arq_texto = open(p_arq_texto, 'w')

	cabecalho = 0
	while (cabecalho < dimensoes + 1):
		numero = arq_binario.read(4)

		if numero == '':
			break;

		numero = struct.unpack('i', numero)[0]

		arq_texto.write(str(numero) + ' ')
		cabecalho += 1

	arq_texto.write('\n')

	cabecalho = 0
	while (cabecalho < dimensoes):
		numero = arq_binario.read(4)

		if not numero:
			break

		numero = struct.unpack('i', numero)[0]

		arq_texto.write(str(numero) + ' ')

		if cabecalho + 1 == dimensoes:
			arq_texto.write('\n')			
			cabecalho = 0
		else:
			cabecalho += 1

	arq_binario.close()
	arq_texto.close()

	print('Fim conversão binario -> texto. Arquivo ' + p_arq_texto + ' criado.')


def converter_novocubo_texto_pra_binario(p_arq_binario, p_arq_texto):
	arq_binario = open(p_arq_binario, 'wb')
	arq_texto = open(p_arq_texto)

	num_linha = 1

	while True:
		linha = arq_texto.readline()

		if linha == '':
			break

		todos_numeros = linha.split(' ')
		num = todos_numeros.pop()

		try:
			int(num)
			print('ERRO na linha ' + str(num_linha) + '!')
		except:
			pass


		for n in todos_numeros:
			value = n.replace('\n', '')
			arq_binario.write(struct.pack('i', int(value)))

		num_linha += 1

	arq_binario.close()
	arq_texto.close()

	print('Fim conversão texto -> binario. Arquivo ' + p_arq_binario + ' criado.')


def converter_cubo():
	p_arq_binario = ''
	p_arq_texto = ''

	path = argv[2]
	dimensoes = int(argv[3])	

	if '.txt' in path:
		p_arq_binario = path.replace('.txt', '.bin')		
		p_arq_texto = path
		converter_novocubo_texto_pra_binario(p_arq_binario, p_arq_texto)
	elif '.bin' in path:
		p_arq_binario = path		
		p_arq_texto = path.replace('.bin', '.txt')
		converter_novocubo_binario_pra_texto(p_arq_binario, p_arq_texto, dimensoes)		
	else:
		raise Exception('Extensao do arquivo nao identificada...')


def ler_binario(path, tid, dimensoes):
	arq = open(path, 'rb')

	t = 1
	tupla = ''

	arq.seek((dimensoes + 1) * 4, 1)

	while t < tid:
		arq.seek(dimensoes * 4, 1)
		t += 1

	for d in range(dimensoes):
		a = arq.read(4)[0]

		if a == '':
			arq.close()
			return

		tupla += str(a) + ' '

	tupla = tupla[:-1]

	print("A tupla\n\n" + str(tid) +":(" + tupla + ")\n\nfoi copiada para a clipboard!")

	pyperclip.copy(tupla)
	arq.close()

	return tupla


def ler_texto(path, tid):
	arq = open(path)

	t = 1

	arq.readline()

	while t <= tid:
		tupla = arq.readline()[:-2]

		if tupla == '':
			arq.close()
			return

		t += 1

	print("A tupla\n\n" + str(tid) +":(" + tupla + ")\n\nfoi copiada para a clipboard!")

	pyperclip.copy(tupla)
	arq.close()

	return tupla


def ler():
	path = argv[2]
	dimensoes = int(argv[3])
	tid = int(argv[4])

	if '.bin' in path:
		return ler_binario(path, tid, dimensoes)
	elif '.txt' == 'path':
		return ler_texto(path, tid)

	return None


def limpar_tudo():
	resposta = input('\nVocê tem certeza que quer limpar todos os cubos? <s/N>')

	if resposta != 's' and resposta != 'S':
		return

	system('del novocubo_*')
	system('dir')
	print("\nTodos os cubos foram limpos!")


def exibir_ajuda():
	print('\nGERAR CUBO:')
	print('\"gerar <cardinalidade> <dimensões> <tuplas> <tid> <[binario] ou [texto]>\"')

	print('\nCONVERTER CUBO:')
	print('\"converter <path_base> <dimensões>\"')	

	print('\nLER CUBO:')
	print('\"ler <path_base> <tid> <dimensões>\"')	

	print('\nDELETAR CUBO:')
	print('\"limpar\"')		


if __name__ == '__main__':
	system('cls')
	try:
		if argv[1] == 'gerar':
			gerar_cubo()
		elif argv[1] == 'converter':
			converter_cubo()
		elif argv[1] == 'ajuda':
			exibir_ajuda()
		elif argv[1] == 'ler':
			ler()
		elif argv[1] == 'limpar':
			limpar_tudo()

	except Exception as e:
		print(e)
		print('\n\nUm erro aconteceu!')
		exibir_ajuda()