import boxsdk as box
import sys
import os
import pickle

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

file_list = os.listdir('/scratch/ravonklar/')
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
                upl_file = client.folder(item.id).upload(f'/scratch/ravonklar/{filename}')
                print(f'Uploaded {upl_file.name} with ID {upl_file.id}')
                break
    else:
        subitems = client.folder(haloids[id]).get_items()
        itemlist = []
        for item in subitems:
            itemlist.append(item.name)
        if filename not in itemlist:
            upl_file = client.folder(haloids[id]).upload(f'/scratch/ravonklar/{filename}')
            print(f'Uploaded {upl_file.name} with ID {upl_file.id}')



