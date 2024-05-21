import subprocess
import os
import sys
import numpy as np
import pickle

import pathlib
path = pathlib.Path().resolve()

rayfilepath = '/scratch/ravonklar/'

'''
halofiles = os.listdir('/home/ravonklar/setup/TNG100-1/cutouts/')
haloids = []
for filepath in halofiles:
	filepath = filepath.replace('cutout_', '')
	filepath = filepath.replace('.hdf5', '')
	haloids.append(filepath)
'''

inactiveHalos = []
with open("inactiveHalos.pkl", "rb") as f:
	inactiveHalos = pickle.load(f)
	haloid = list(inactiveHalos.keys())[0]
	startInd = inactiveHalos[haloid]
	inactiveHalos.pop(haloid)
with open("inactiveHalos.pkl", "wb") as f:
	pickle.dump(inactiveHalos, f)

print(f"Working on halo {haloid}")

inc = 100

finalGoal = (int(startInd/inc)+1)*inc

'''if startInd >= 4999:
	finalGoal = 10000
else:
	finalGoal = 5000
	'''
inc2 = 10
estimatedNumIncrements = int((finalGoal-startInd-1)/inc2)+1

startingGoal = (int(startInd/inc2)+1)*inc2

print(f"Starting at index {startInd} towards {finalGoal}")

for i in range(estimatedNumIncrements):
	intGoal = startingGoal+inc2*i

	proc = subprocess.Popen(f"cd {path} && python raycontroller.py {haloid} {startInd} {intGoal-startInd}", shell=True, universal_newlines=True)
	proc.wait()
	startInd = intGoal
print(f'INCREMENT {int(finalGoal/inc)} SUCCESSFULLY COMPLETED')

if finalGoal != 10000:
	with open("inactiveHalos.pkl", "wb") as f:
		intdict = {}
		intdict[haloid] = finalGoal
		inactiveHalos = {**intdict, **inactiveHalos}
		pickle.dump(inactiveHalos, f)
else:
	doneHalos = []
	with open("doneHalos.pkl", "rb") as f:
		doneHalos = pickle.load(f)
		doneHalos.append(haloid)
	with open("doneHalos.pkl", "wb") as f:
		pickle.dump(doneHalos, f)
	import boxsdk as box
	import discord
	from discord.ext import commands
	import asyncio
	from dotenv import load_dotenv

	load_dotenv()
	client_id = os.getenv('clientID')
	client_secret = os.getenv('clientSecret')

	oauth = box.OAuth2(
		client_id=client_id,
		client_secret=client_secret
	)

	auth_url, csrf_token = oauth.get_authorization_url(None)

	print(auth_url)

	print('WAITING FOR AUTH CODE')

	token = os.getenv('discordToken')

	auth_code = []

	async def mainBot():
		intents = discord.Intents.default()
		intents.message_content = True
		bot = commands.Bot(command_prefix = ['!'], intents=intents)
		@bot.event
		async def on_ready():
			me = await bot.fetch_user(198200248634572811)
			await me.send(auth_url)
		
		@bot.command(name='code')
		async def receive_code(ctx, code):
			auth_code.append(code)
			me = await bot.fetch_user(198200248634572811)
			await me.send('Received!')
			await bot.close()
		
		return bot

	bot = asyncio.run(mainBot())
	bot.run(token)

	access_token, refresh_token = oauth.authenticate(auth_code[0])
	client = box.Client(oauth)

	user = client.user().get()
	print(f'User ID is {user.id}')

	folder_id = 237212492701
	items = client.folder(folder_id).get_items()

	print(items)

	file_list = os.listdir(rayfilepath)
	haloids = {}
	for item in items:
		haloids[item.name] = item.id

	for filename in file_list:
		idname = filename.replace('gal', '')
		id = idname.split('_',1)[0]
		if id not in haloids:
			print(f"NEW ID: {id}")
			print(haloids)
			for item in items:
				print(item.name)
				if item.name == id:
					haloids[id] = item.id
					print("FOUND")
					break
		if id not in haloids:
			print(f"NOT FOUND: {id}")
			print(haloids)
			client.folder(folder_id).create_subfolder(id)
			items = client.folder(folder_id).get_items()
			for item in items:
				if item.name == id:
					haloids[id] = item.id
					upl_file = client.folder(item.id).upload(f'{rayfilepath}{filename}')
					print(f'Uploaded {upl_file.name} with ID {upl_file.id}')
					break
		else:
			subitems = client.folder(haloids[id]).get_items()
			itemlist = []
			for item in subitems:
				itemlist.append(item.name)
			if filename not in itemlist:
				upl_file = client.folder(haloids[id]).upload(f'{rayfilepath}{filename}')
				print(f'Uploaded {upl_file.name} with ID {upl_file.id}')
	'''
	if (len(doneHalos)-1)%6 == 0:
		#start rayuploader
		proc = subprocess.Popen(f"cd {path} && sbatch uploadscript.bat")
		proc.wait()
		print('UPLOADING')
		'''

'''
if len(inactiveHalos) > 0:
	proc = subprocess.Popen(f"cd {path} && sbatch bscrip.bat")
	proc.wait()
'''